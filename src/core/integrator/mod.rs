//! Integrator

#![allow(dead_code)]

mod common;
mod sampler_integrator;

use crate::core::scene::Scene;
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
}

/// Atomic reference counted `Integrator`.
pub type ArcIntegrator = Arc<dyn Integrator + Send + Sync>;
