//! Grid Density Media

use bumpalo::Bump;
use core::geometry::*;
use core::interaction::MediumInteraction;
use core::medium::HenyeyGreenstein;
use core::medium::Medium;
use core::pbrt::*;
use core::sampler::*;
use core::spectrum::*;
use std::sync::Arc;

/// Implements medium densities at a regular 3D grid of positions, similar to the
/// way that the ImageTexture represents images with a 2D grid of samples. These
/// samples are interpolated to compute the density at positions between the sample
/// points.
pub struct GridDensityMedium {
    /// Absorption cross section `σa` is the probability density that light is
    /// absorbed per unit distance traveled in the medium
    _sigma_a: Spectrum,

    /// Scattering coefficient `σs` is the probability of an out-scattering
    /// event occurring per unit distance
    sigma_s: Spectrum,

    /// Total reduction in radiance due to absorption and out-scattering
    /// `σt = σs + σa`. This combined effect of absorption and out-scattering is
    /// called attenuation or extinction.
    sigma_t: Float,

    /// The asymmetry parameter for Henyey-Greenstein phase function.
    g: Float,

    /// Grid size in x-direction.
    nx: usize,

    /// Grid size in y-direction.
    ny: usize,

    /// Grid size in z-direction.
    nz: usize,

    /// Transformaton from world-space to medium-space.
    world_to_medium: Transform,

    /// Density values in the grid.
    density: Vec<Float>,

    /// 1 / maximum density value in the grid.
    inv_max_density: Float,
}

impl GridDensityMedium {
    /// Create a new `GridDensityMedium `.
    ///
    /// * `sigma_a`         - Absorption cross section `σa`.
    /// * `sigma_s`         - Scattering coefficient `σs`.
    /// * `g`               - The asymmetry parameter for Henyey-Greenstein phase
    ///                       function.
    /// * `nx`              - Grid size in x-direction.
    /// * `ny`              - Grid size in y-direction.
    /// * `nz`              - Grid size in z-direction.
    /// * `medium_to_world` - Medium-space to world-space transformation.
    /// * `density`         - Density values in the grid.
    pub fn new(
        sigma_a: Spectrum,
        sigma_s: Spectrum,
        g: Float,
        nx: usize,
        ny: usize,
        nz: usize,
        medium_to_world: Transform,
        density: Vec<Float>,
    ) -> Self {
        // Precompute values for Monte Carlo sampling of `GridDensityMedium`.
        let sigma_t = sigma_s + sigma_a;
        if Spectrum::new(sigma_t[0]) != sigma_t {
            error!("GridDensityMedium requires a spectrally uniform attenuation coefficient!");
        }

        let max_density = density.iter().fold(0.0, |a, &x| max(a, x));
        let inv_max_density = 1.0 / max_density;
        Self {
            _sigma_a: sigma_a,
            sigma_s,
            sigma_t: sigma_t[0],
            g,
            nx,
            ny,
            nz,
            world_to_medium: medium_to_world.inverse(),
            density,
            inv_max_density,
        }
    }

    /// Reconstruct the volume density function at the given position.
    ///
    /// * `p` - Sample position.
    fn density(&self, p: &Point3f) -> Float {
        // Compute voxel coordinates and offsets for `p`.
        let p_samples = Point3f::new(
            p.x * self.nx as Float - 0.5,
            p.y * self.ny as Float - 0.5,
            p.z * self.nz as Float - 0.5,
        );
        let pi = Point3i::from(p_samples.floor());
        let d = p_samples - Point3f::from(pi);

        // Trilinearly interpolate density values to compute local density.
        let d00 = lerp(d.x, self.d(&pi), self.d(&(pi + Vector3i::new(1, 0, 0))));
        let d10 = lerp(
            d.x,
            self.d(&(pi + Vector3i::new(0, 1, 0))),
            self.d(&(pi + Vector3i::new(1, 1, 0))),
        );
        let d01 = lerp(
            d.x,
            self.d(&(pi + Vector3i::new(0, 0, 1))),
            self.d(&(pi + Vector3i::new(1, 0, 1))),
        );
        let d11 = lerp(
            d.x,
            self.d(&(pi + Vector3i::new(0, 1, 1))),
            self.d(&(pi + Vector3i::new(1, 1, 1))),
        );
        let d0 = lerp(d.y, d00, d10);
        let d1 = lerp(d.y, d01, d11);
        lerp(d.z, d0, d1)
    }

    /// Returns the density at the given integer sample position.
    ///
    /// * `p` - Sample position.
    fn d(&self, p: &Point3i) -> Float {
        let sample_bounds = Bounds3i::new(
            Point3i::ZERO,
            Point3i::new(self.nx as Int, self.ny as Int, self.nz as Int),
        );
        if !sample_bounds.contains_exclusive(&p) {
            0.0
        } else {
            let i = (p.z * self.ny as Int + p.y) * self.nx as Int + p.x;
            self.density[i as usize]
        }
    }
}

impl Medium for GridDensityMedium {
    /// Returns the beam transmittance along a given ray.
    ///
    /// * `ray`     - The ray.
    /// * `sampler` - The sampler.
    fn tr(&self, ray: &Ray, sampler: &mut ArcSampler) -> Spectrum {
        let r = Ray::new(
            ray.o,
            ray.d.normalize(),
            ray.t_max * ray.d.length(),
            0.0,
            None,
        );
        let ray_medium = self.world_to_medium.transform_ray(&r);

        // Compute [t_min, t_max] interval of `ray`'s overlap with medium bounds.
        let b = Bounds3f::new(Point3f::ZERO, Point3f::new(1.0, 1.0, 1.0));
        if let Some((t_min, t_max)) = b.intersect_p(&ray_medium) {
            let sampler = Arc::get_mut(sampler).unwrap();

            // Perform ratio tracking to estimate the transmittance value.
            let mut tr = 1.0;
            let mut t = t_min;
            loop {
                t -= (1.0 - sampler.get_1d()).ln() * self.inv_max_density / self.sigma_t;
                if t >= t_max {
                    break;
                }

                let density = self.density(&ray_medium.at(t));
                tr *= 1.0 - max(0.0, density * self.inv_max_density);

                // Added after book publication: when transmittance gets low,
                // start applying Russian roulette to terminate sampling.
                const RR_THRESHOLD: Float = 0.1;
                if tr < RR_THRESHOLD {
                    let q = max(0.05, 1.0 - tr);
                    if sampler.get_1d() < q {
                        return Spectrum::ZERO;
                    }
                    tr /= 1.0 - q;
                }
            }
            Spectrum::new(tr)
        } else {
            Spectrum::ONE
        }
    }

    /// Samples a medium scattering interaction along a world-space ray.
    ///
    /// The ray will generally have been intersected against the scene geometry;
    /// thus, implementations of this method shouldn’t ever sample a medium
    /// interaction at a point on the ray beyond its `t_max` value.
    ///
    /// NOTE: Calling code will need to assign this medium as we cannot pass
    /// back and `ArcMedium` out of here for `Self`.
    ///
    /// * `arena`   - The arena for memory allocations.
    /// * `ray`     - The ray.
    /// * `sampler` - The sampler.
    fn sample<'arena>(
        &self,
        arena: &'arena Bump,
        ray: &Ray,
        sampler: &mut ArcSampler,
    ) -> (Spectrum, Option<MediumInteraction<'arena>>) {
        let r = Ray::new(
            ray.o,
            ray.d.normalize(),
            ray.t_max * ray.d.length(),
            0.0,
            None,
        );
        let ray_medium = self.world_to_medium.transform_ray(&r);

        // Compute [t_min, t_max] interval of `ray`'s overlap with medium bounds.
        let b = Bounds3f::new(Point3f::ZERO, Point3f::new(1.0, 1.0, 1.0));
        if let Some((t_min, t_max)) = b.intersect_p(&ray_medium) {
            // Run delta-tracking iterations to sample a medium interaction.
            let sampler = Arc::get_mut(sampler).unwrap();
            let mut t = t_min;
            loop {
                t -= (1.0 - sampler.get_1d()).ln() * self.inv_max_density / self.sigma_t;
                if t >= t_max {
                    break;
                }

                if self.density(&ray_medium.at(t)) * self.inv_max_density > sampler.get_1d() {
                    // Populate `mi` with medium interaction information and return.
                    let phase = HenyeyGreenstein::alloc(arena, self.g);
                    let mi = MediumInteraction::new(ray.at(t), -ray.d, ray.time, None, phase);
                    return (self.sigma_s / self.sigma_t, Some(mi));
                }
            }
            (Spectrum::ONE, None)
        } else {
            (Spectrum::ONE, None)
        }
    }
}
