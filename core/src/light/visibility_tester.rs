//! Visibility Tester

use crate::interaction::*;
use crate::sampler::*;
use crate::scene::*;
use crate::spectrum::*;

/// VisibilityTester allows lights to return a radiance value under the
/// assumption that the reference point and light source are mutually
/// visible.
#[derive(Clone)]
pub struct VisibilityTester {
    /// One endpoint of shadow ray.
    pub p0: Hit,

    /// Second endpoint of shadow ray.
    pub p1: Hit,
}

impl VisibilityTester {
    /// Create a new `VisibilityTester` for given endpoints of a shadow ray.
    ///
    /// * `p0` - One endpoint of shadow ray.
    /// * `p1` - Second endpoint of shadow ray.
    pub fn new(p0: Hit, p1: Hit) -> Self {
        Self { p0, p1 }
    }

    /// Traces a shadow ray between `p0` and `p1` through the scene and returns
    /// true if the points are visible to each other.
    ///
    /// * `scene` - The scene.
    pub fn unoccluded(&self, scene: &Scene) -> bool {
        !scene.intersect_p(&self.p0.spawn_ray_to_hit(&self.p1))
    }

    /// Computes the beam transmittance, the fraction of radiance transmitted
    /// along the segment between the two points. It accounts for both attenuation
    /// in participating media as well as any surfaces that block the ray completely.
    ///
    /// * `scene`   - The scene.
    /// * `sampler` - The sampler.
    pub fn tr(&self, scene: &Scene, sampler: &mut ArcSampler) -> Spectrum {
        let mut ray = self.p0.spawn_ray_to_hit(&self.p1);
        let mut tr_val = Spectrum::ONE;

        loop {
            let isect = scene.intersect(&mut ray);

            // Handle opaque surface along ray's path.
            if let Some(isect) = isect.as_ref() {
                if let Some(primitive) = isect.primitive {
                    if primitive.get_material().is_some() {
                        return Spectrum::ZERO;
                    }
                }
            }

            // Update transmittance for current ray segment.
            if let Some(medium) = ray.medium.as_ref() {
                tr_val *= medium.tr(&ray, sampler);
            }
            // Generate next ray segment or return final transmittance.
            if let Some(isect) = isect {
                ray = isect.hit.spawn_ray_to_hit(&self.p1);
            } else {
                break;
            }
        }

        tr_val
    }
}
