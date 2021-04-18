//! Visibility Tester

use crate::core::geometry::*;
use crate::core::sampler::*;
use crate::core::scene::*;
use crate::core::spectrum::*;
use std::sync::Arc;

/// VisibilityTester allows lights to return a radiance value under the
/// assumption that the reference point and light source are mutually
/// visible.
#[derive(Clone)]
pub struct VisibilityTester {
    /// One endpoint of shadow ray.
    pub p0: ArcInteraction,

    /// Second endpoint of shadow ray.
    pub p1: ArcInteraction,
}

impl VisibilityTester {
    /// Create a new `VisibilityTester` for given endpoints of a shadow ray.
    ///
    /// * `p0` - One endpoint of shadow ray.
    /// * `p1` - Second endpoint of shadow ray.
    pub fn new(p0: ArcInteraction, p1: ArcInteraction) -> Self {
        Self {
            p0: p0.clone(),
            p1: p1.clone(),
        }
    }

    /// Traces a shadow ray between `p0` and `p1` through the scene and returns
    /// true if the points are visible to each other.
    ///
    /// * `scene` - The scene.
    pub fn unoccluded(&self, scene: Arc<Scene>) -> bool {
        scene.intersect_p(&self.p0.get_hit().spawn_ray_to(&self.p1.get_hit().p))
    }

    /// Computes the beam transmittance, the fraction of radiance transmitted
    /// along the segment between the two points. It accounts for both attenuation
    /// in participating media as well as any surfaces that block the ray completely.
    ///
    /// * `scene`   - The scene.
    /// * `sampler` - The sampler.
    pub fn tr(&self, scene: Arc<Scene>, sampler: ArcSampler) -> Spectrum {
        let mut ray = self.p0.get_hit().spawn_ray_to(&self.p1.get_hit().p);
        let mut tr = Spectrum::new(1.0);

        loop {
            if let Some(isect) = scene.intersect(&mut ray) {
                // Handle opaque surface along ray's path.
                if let Some(_material) = isect.primitive.map(|p| p.get_material()) {
                    return Spectrum::new(0.0);
                }

                // Update transmittance for current ray segment.
                let medium = ray.medium.clone();
                if let Some(tr2) = medium.map(|medium| medium.tr(&ray, sampler.clone())) {
                    tr *= tr2;
                }

                // Generate next ray segment or return final transmittance.
                ray = isect.get_hit().spawn_ray_to(&self.p1.get_hit().p);
            } else {
                // Update transmittance for current ray segment.
                let medium = ray.medium.clone();
                if let Some(tr2) = medium.map(|medium| medium.tr(&ray, sampler.clone())) {
                    tr *= tr2;
                }
                break;
            }
        }

        tr
    }
}
