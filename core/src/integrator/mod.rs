//! Integrator

mod common;
mod sampler_integrator;

use crate::geometry::*;
use crate::sampler::*;
use crate::scene::Scene;
use crate::spectrum::*;

// Re-export.
pub use common::*;
pub use sampler_integrator::*;

/// Integrator interface.
pub trait Integrator {
    /// Render the scene.
    ///
    /// * `scene` - The scene.
    fn render(&self, scene: &Scene);

    /// Preprocess the scene.
    ///
    /// * `scene` - The scene
    fn preprocess(&mut self, scene: &Scene);

    /// Returns the incident radiance at the origin of a given ray.
    ///
    /// * `ray`     - The ray.
    /// * `scene`   - The scene.
    /// * `sampler` - The sampler.
    /// * `depth`   - The recursion depth.
    fn li(&self, _ray: &mut Ray, _scene: &Scene, _sampler: &mut dyn Sampler, _depth: usize) -> Spectrum {
        Spectrum::ZERO
    }

    /// Returns the cropped pixel bounds of the image.
    fn get_cropped_pixel_bounds(&self) -> Bounds2i;
}
