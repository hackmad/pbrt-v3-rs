//! Integrator

mod sampler_integrator;
mod common;

use crate::geometry::*;
use crate::sampler::*;
use crate::scene::Scene;
use crate::spectrum::*;
use std::sync::Arc;

// Re-export.
pub use common::*;
pub use sampler_integrator::*;

/// Integrator interface.
pub trait Integrator {
    /// Render the scene.
    ///
    /// * `scene` - The scene.
    fn render(&mut self, scene: Arc<Scene>);

    /// Returns the incident radiance at the origin of a given ray.
    ///
    /// * `ray`     - The ray.
    /// * `scene`   - The scene.
    /// * `sampler` - The sampler.
    /// * `depth`   - The recursion depth.
    fn li(
        &self,
        _ray: &mut Ray,
        _scene: Arc<Scene>,
        _sampler: &mut ArcSampler,
        _depth: usize,
    ) -> Spectrum {
        Spectrum::new(0.0)
    }
}

/// Atomic reference counted `Integrator`.
pub type ArcIntegrator = Arc<dyn Integrator + Send + Sync>;
