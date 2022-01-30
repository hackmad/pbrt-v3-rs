//! Integrator

mod common;
mod sampler_integrator;

use crate::geometry::*;
use crate::sampler::*;
use crate::scene::Scene;
use crate::spectrum::*;
use bumpalo::Bump;

// Re-export.
pub use common::*;
pub use sampler_integrator::*;

/// Integrator interface.
pub trait Integrator {
    /// Render the scene.
    ///
    /// * `scene` - The scene.
    fn render(&self, scene: &Scene);

    /// Returns the incident radiance at the origin of a given ray.
    ///
    /// * `arena`   - The arena for memory allocations.
    /// * `ray`     - The ray.
    /// * `scene`   - The scene.
    /// * `sampler` - The sampler.
    /// * `depth`   - The recursion depth.
    fn li(
        &self,
        _arena: &Bump,
        _ray: &mut Ray,
        _scene: &Scene,
        _sampler: &mut ArcSampler,
        _depth: usize,
    ) -> Spectrum {
        Spectrum::ZERO
    }
}
