//! Textures

#![allow(dead_code)]
use crate::geometry::*;
use crate::pbrt::Float;
use crate::spectrum::Spectrum;
use std::collections::HashMap;
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

/// Map of floating point textures.
pub type FloatTextureMap = HashMap<String, ArcTexture<Float>>;

/// Map of spectrum textures.
pub type SpectrumTextureMap = HashMap<String, ArcTexture<Spectrum>>;

mod common;
mod mapping;

// Re-export
pub use common::*;
pub use mapping::*;
