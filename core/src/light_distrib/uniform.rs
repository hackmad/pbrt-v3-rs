//! Uniform Light Distribution.

use super::LightDistribution;
use crate::geometry::*;
use crate::sampling::*;
use crate::scene::*;
use std::sync::Arc;

/// The simplest possible implementation of `LightDistribution`: this returns a uniform distribution over all light
/// sources, ignoring the provided point. This approach works well for very simple scenes, but is quite ineffective for
/// scenes with more than a handful of light sources. (This was the sampling method originally used for the
/// `PathIntegrator` and the `VolPathIntegrator` in the printed book, though without the `UniformLightDistribution`
/// class.)
pub struct UniformLightDistribution {
    distrib: Option<Arc<Distribution1D>>,
}

impl UniformLightDistribution {
    /// Create a new instance of `PowerLightDistribution`.
    ///
    /// * `scene` - The scene.
    pub fn new(scene: &Scene) -> Self {
        let prob = vec![1.0; scene.lights.len()];
        Self {
            distrib: Some(Arc::new(Distribution1D::new(prob))),
        }
    }
}

impl LightDistribution for UniformLightDistribution {
    /// Given a point |p| in space, this method returns a (hopefully effective) sampling distribution for light sources
    /// at that point.
    fn lookup(&self, _p: &Point3f) -> Option<Arc<Distribution1D>> {
        self.distrib.as_ref().map(Arc::clone)
    }
}
