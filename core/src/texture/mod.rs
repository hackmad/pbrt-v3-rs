//! Textures

#![allow(dead_code)]
use crate::geometry::*;
use crate::interaction::*;
use crate::pbrt::Float;
use crate::spectrum::Spectrum;
use std::collections::HashMap;
use std::sync::Arc;

/// Texture interface.
pub trait Texture<T: Copy> {
    /// Evaluate the texture at surface interaction.
    ///
    /// * `hit` - Surface interaction hit.
    /// * `uv`  - Surface interaction uv.
    /// * `der` - Surface interaction derivatives.
    fn evaluate(&self, hit: &Hit, uv: &Point2f, der: &Derivatives) -> T;
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
