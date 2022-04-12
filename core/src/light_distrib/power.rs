//! Power Light Distribution.

use super::LightDistribution;
use crate::geometry::*;
use crate::integrator::compute_light_power_distribution;
use crate::sampling::*;
use crate::scene::*;
use std::sync::Arc;

/// PowerLightDistribution returns a distribution with sampling probability
/// proportional to the total emitted power for each light. (It also ignores
/// the provided point |p|.)  This approach works well for scenes where
/// there the most powerful lights are also the most important contributors
/// to lighting in the scene, but doesn't do well if there are many lights
/// and if different lights are relatively important in some areas of the
/// scene and unimportant in others. (This was the default sampling method
/// used for the BDPT integrator and MLT integrator in the printed book,
/// though also without the PowerLightDistribution class.)
pub struct PowerLightDistribution {
    distrib: Option<Arc<Distribution1D>>,
}

impl PowerLightDistribution {
    /// Create a new instance of `PowerLightDistribution`.
    ///
    /// * `scene` - The scene.
    pub fn new(scene: &Scene) -> Self {
        Self {
            distrib: compute_light_power_distribution(scene).map(Arc::new),
        }
    }
}

impl LightDistribution for PowerLightDistribution {
    /// Given a point |p| in space, this method returns a (hopefully effective)
    /// sampling distribution for light sources at that point.
    fn lookup(&self, _p: &Point3f) -> Option<Arc<Distribution1D>> {
        self.distrib.as_ref().map(Arc::clone)
    }
}
