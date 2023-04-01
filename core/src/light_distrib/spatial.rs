//! Spatial Light Distribution.

use super::LightDistribution;
use crate::geometry::*;
use crate::interaction::Hit;
use crate::light::*;
use crate::low_discrepency::radical_inverse;
use crate::pbrt::*;
use crate::sampling::*;
use crate::scene::*;
use crate::spectrum::*;
use arc_swap::ArcSwapOption;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

/// Voxel coordinates are packed into a usize (64-bit) for hash table lookups; 10 bits are allocated
/// to each coordinate. INVALID_PACKED_POS is an impossible packed coordinate value, which we use to
/// represent.
const INVALID_PACKED_POS: usize = 0xffffffffffffffff;

struct HashEntry {
    packed_pos: AtomicUsize,
    distribution: ArcSwapOption<Distribution1D>,
}

impl Default for HashEntry {
    /// Returns the "default value" for `HashEntry`.
    fn default() -> Self {
        Self {
            packed_pos: AtomicUsize::new(INVALID_PACKED_POS),
            distribution: ArcSwapOption::const_empty(),
        }
    }
}

/// A spatially-varying light distribution that adjusts the probability of sampling a light source
/// based on an estimate of its contribution to a region of space.  A fixed voxel grid is imposed
/// over the scene bounds and a sampling distribution is computed as needed for each voxel.
pub struct SpatialLightDistribution {
    lights: Vec<ArcLight>,
    world_bound: Bounds3f,
    n_voxels: [usize; 3],
    hash_table: Vec<HashEntry>,
}

impl SpatialLightDistribution {
    /// Create a new instance of `SpatialLightDistribution`.
    ///
    /// * `scene`      - The scene.
    /// * `max_voxels` - Maximum number of voxels (default to 64).
    pub fn new(scene: &Scene, max_voxels: usize) -> Self {
        // Compute the number of voxels so that the widest scene bounding box dimension has maxVoxels
        // voxels and the other dimensions have a number of voxels so that voxels are roughly cube shaped.
        let b = scene.world_bound;
        let diag = b.diagonal();
        let bmax = diag[b.maximum_extent()];
        let mut n_voxels = [0; 3];
        for i in 0..3 {
            n_voxels[i] = max(1_usize, (diag[i] / bmax * max_voxels as Float).round() as usize) as usize;

            // In the Lookup() method, we require that 20 or fewer bits be sufficient to represent
            // each coordinate value. It's fairly hard to imagine that this would ever be a problem.
            assert!(n_voxels[i] < (1 << 20));
        }

        info!(
            "SpatialLightDistribution: scene bounds {}, voxel res ({}, {}, {})",
            b, n_voxels[0], n_voxels[1], n_voxels[2]
        );

        let hash_table_size = 4 * n_voxels[0] * n_voxels[1] * n_voxels[2];
        Self {
            lights: scene.lights.iter().map(Arc::clone).collect(),
            world_bound: b,
            n_voxels,
            hash_table: (0..hash_table_size).map(|_| HashEntry::default()).collect(),
        }
    }

    // Compute the sampling distribution for the voxel with integer coordiantes given by "pi".
    fn compute_distribution(&self, pi: &Point3i) -> Distribution1D {
        // Compute the world-space bounding box of the voxel corresponding to |pi|.
        let p0 = Point3f::new(
            pi[0] as Float / self.n_voxels[0] as Float,
            pi[1] as Float / self.n_voxels[1] as Float,
            pi[2] as Float / self.n_voxels[2] as Float,
        );
        let p1 = Point3f::new(
            (pi[0] + 1) as Float / self.n_voxels[0] as Float,
            (pi[1] + 1) as Float / self.n_voxels[1] as Float,
            (pi[2] + 1) as Float / self.n_voxels[2] as Float,
        );
        let voxel_bounds = Bounds3f::new(self.world_bound.lerp(&p0), self.world_bound.lerp(&p1));

        // Compute the sampling distribution. Sample a number of points inside voxelBounds using a
        // 3D Halton sequence; at each one, sample each light source and compute a weight based on
        // Li/pdf for the light's sample (ignoring visibility between the point in the voxel and the
        // point on the light source) as an approximation to how much the light is likely to contribute
        // to illumination in the voxel.
        const N_SAMPLES: usize = 128;
        let n_lights = self.lights.len();
        let mut light_contrib = vec![Float::default(); n_lights];
        for i in 0..N_SAMPLES {
            let po = voxel_bounds.lerp(&Point3f::new(
                radical_inverse(0_u16, i as u64),
                radical_inverse(1_u16, i as u64),
                radical_inverse(2_u16, i as u64),
            ));
            let intr = Hit::new(
                po,
                0.0, /* time */
                Vector3f::ZERO,
                Vector3f::new(1.0, 0.0, 0.0),
                Normal3f::ZERO,
                None,
            );

            // Use the next two Halton dimensions to sample a point on the
            // light source.
            let u = Point2f::new(radical_inverse(3_u16, i as u64), radical_inverse(4_u16, i as u64));
            for (j, light) in self.lights.iter().enumerate() {
                if let Some(li) = light.sample_li(&intr, &u) {
                    if li.pdf > 0.0 {
                        // TODO: Look at tracing shadow rays / computing beam transmittance. Probably
                        // shouldn't give those full weight but instead e.g. have an occluded shadow
                        // ray scale down the contribution by 10 or something.
                        light_contrib[j] += li.value.y() / li.pdf;
                    }
                }
            }
        }

        // We don't want to leave any lights with a zero probability; it's possible that a light
        // contributes to points in the voxel even though we didn't find such a point when sampling
        // above. Therefore, compute a minimum (small) weight and ensure that all lights are given at
        // least the corresponding probability.
        let sum_contrib: Float = light_contrib.iter().sum();
        let avg_contrib = sum_contrib / (N_SAMPLES * light_contrib.len()) as Float;
        let min_contrib = if avg_contrib > 0.0 { 0.001 * avg_contrib } else { 1.0 };
        for (i, contrib) in light_contrib.iter_mut().enumerate() {
            debug!("Voxel pi = {pi}, light {i} contrib = {contrib}");
            *contrib = max(*contrib, min_contrib);
        }
        info!("Initialized light distribution in voxel pi = {pi}, avgContrib = {avg_contrib}",);

        // Compute a sampling distribution from the accumulated contributions.
        Distribution1D::new(light_contrib)
    }
}

impl LightDistribution for SpatialLightDistribution {
    /// Given a point |p| in space, this method returns a (hopefully effective) sampling distribution
    /// for light sources at that point.
    fn lookup(&self, p: &Point3f) -> Option<Arc<Distribution1D>> {
        // First, compute integer voxel coordinates for the given point |p| with respect to the overall
        // voxel grid.
        let offset = self.world_bound.offset(p); // offset in [0,1].
        let mut pi = Point3i::ZERO;
        for i in 0..3 {
            // The clamp should almost never be necessary, but is there to be robust to computed
            // intersection points being slightly outside the scene bounds due to floating-point
            // roundoff error.
            pi[i] = clamp(
                (offset[i] * self.n_voxels[i] as Float) as Int,
                0,
                self.n_voxels[i] as Int - 1,
            );
        }

        let hash_table_size = self.hash_table.len();

        // Pack the 3D integer voxel coordinates into a single 64-bit value.
        let packed_pos = ((pi[0] as usize) << 40) | ((pi[1] as usize) << 20) | pi[2] as usize;
        assert_ne!(packed_pos, INVALID_PACKED_POS);

        // Compute a hash value from the packed voxel coordinates. We could just take packedPos mod
        // the hash table size, but since packedPos isn't necessarily well distributed on its own,
        // it's worthwhile to do a little work to make sure that its bits values are individually
        // fairly random. For details of and motivation for the following, see:
        // http://zimbry.blogspot.ch/2011/09/better-bit-mixing-improving-on.html
        let mut hash = packed_pos;
        hash ^= hash >> 31;
        hash *= hash.wrapping_mul(0x7fb5d329728ea185);
        hash ^= hash >> 27;
        hash *= hash.wrapping_mul(0x81dadef4bc2dd44d);
        hash ^= hash >> 33;
        hash %= hash_table_size;

        // Now, see if the hash table already has an entry for the voxel. We'll use quadratic probing
        // when the hash table entry is already used for another value; step stores the square root
        // of the probe step.
        let mut step = 1;
        loop {
            let entry = &self.hash_table[hash];

            // Does the hash table entry at offset |hash| match the current point?
            let entry_packed_pos = entry.packed_pos.load(Ordering::Acquire);
            if entry_packed_pos == packed_pos {
                let dist = entry.distribution.load().as_ref().map(Arc::clone);
                break dist;
            } else if entry_packed_pos != INVALID_PACKED_POS {
                // The hash table entry we're checking has already been allocated for another voxel.
                // Advance to the next entry with quadratic probing.
                hash += step * step;
                if hash >= hash_table_size {
                    hash %= hash_table_size;
                }
                step += 1;
            } else {
                // We have found an invalid entry. (Though this may have changed since the load into
                // entryPackedPos above.) Use an atomic compare/exchange to try to claim this entry
                // for the current position.
                let invalid = INVALID_PACKED_POS;
                if entry
                    .packed_pos
                    .compare_exchange_weak(invalid, packed_pos, Ordering::Acquire, Ordering::Relaxed)
                    .is_ok()
                {
                    // Success; we've claimed this position for this voxel's distribution. Now compute
                    // the sampling distribution and add it to the hash table. As long as packedPos
                    // has been set but the entry's distribution pointer is nullptr, any other threads
                    // looking up the distribution for this voxel will spin wait until the distribution
                    // pointer is written.
                    let dist = Arc::new(self.compute_distribution(&pi));
                    entry.distribution.swap(Some(Arc::clone(&dist)));
                    break Some(dist);
                }
            }
        }
    }
}
