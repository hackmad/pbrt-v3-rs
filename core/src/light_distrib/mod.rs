//! Light Distribution.

mod power;
mod spatial;
mod uniform;

pub use power::*;
pub use spatial::*;
pub use uniform::*;

use crate::geometry::*;
use crate::sampling::*;
use crate::scene::*;
use std::sync::Arc;

/// Light sampling strategy.
#[derive(Copy, Clone)]
pub enum LightSampleStategy {
    /// Sample all light sources uniformly.
    Uniform,

    /// Samples light sources according to their emitted power.
    Power,

    /// Compute light contributions in regions of the scene and samples from a
    /// related distribution.
    Spatial,
}

impl From<&str> for LightSampleStategy {
    /// Returns a `LightSampleStrategy` given a string name.
    fn from(name: &str) -> Self {
        match name {
            "uniform" => Self::Uniform,
            "power" => Self::Power,
            "spatial" => Self::Spatial,
            _ => {
                error!(
                    "Light sample distribution type '{}' unknown. Using 'spatial'.",
                    name
                );
                Self::Spatial
            }
        }
    }
}

/// Interface of light distribution implementations that provide probability
/// distributions for sampling light sources at a given point in space.
pub trait LightDistribution {
    /// Given a point |p| in space, this method returns a (hopefully effective)
    /// sampling distribution for light sources at that point.
    fn lookup(&self, p: &Point3f) -> Option<Arc<Distribution1D>>;
}

/// Atomic reference counted `LightDistribution `.
pub type ArcLightDistribution = Arc<dyn LightDistribution + Send + Sync>;

/// Returns a smart pointer to a new `LightDistribution` implementation.
///
/// * `strategy` - The strategy to use for light sampling.
/// * `scene`    - The scene.
pub fn create_light_sample_distribution(
    strategy: LightSampleStategy,
    scene: &Scene,
) -> ArcLightDistribution {
    let strategy = if scene.lights.len() == 1 {
        LightSampleStategy::Uniform
    } else {
        strategy
    };
    match strategy {
        LightSampleStategy::Uniform => Arc::new(UniformLightDistribution::new(scene)),
        LightSampleStategy::Power => Arc::new(PowerLightDistribution::new(scene)),
        LightSampleStategy::Spatial => Arc::new(SpatialLightDistribution::new(scene, 64)),
    }
}
