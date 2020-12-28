//! Textures

#![allow(dead_code)]
use crate::core::geometry::*;
use std::sync::Arc;

/// Texture interface.
pub trait Texture<T: Copy> {
    /// Evaluate the texture at surface interaction.
    ///
    /// * `si` - Surface interaction.
    fn evaluate(&self, si: &SurfaceInteraction) -> T;
}

/// Atomic reference counted `Texture`.
pub type ArcTexture<T> = Arc<dyn Texture<T> + Send + Sync>;

mod common;
mod mapping;

// Re-export
pub use common::*;
pub use mapping::*;
