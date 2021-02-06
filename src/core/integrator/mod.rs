//! Integrator

use crate::core::scene::*;
use std::sync::Arc;

/// Integrator interface.
pub trait Integrator {
    /// Render the scene.
    ///
    /// * `scene` - The scene.
    fn render(&self, scene: Arc<Scene>);
}

/// Atomic reference counted `Integrator`.
pub type ArcIntegrator = Arc<dyn Integrator + Send + Sync>;
