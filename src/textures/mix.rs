//! Mix Texture

#[allow(dead_code)]
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::texture::*;
use std::ops::{Add, Mul};

/// Implements a texture that linearly interpolates between two textures with a
/// third texture.
#[derive(Clone)]
pub struct MixTexture<T> {
    /// First texture.
    tex1: ArcTexture<T>,

    /// Second texture.
    tex2: ArcTexture<T>,

    /// Scale amount.
    amount: ArcTexture<Float>,
}

impl<T> MixTexture<T> {
    /// Create a new `MixTexture<T>`.
    ///
    /// * `tex1`   - The first texture.
    /// * `tex2`   - The second texture.
    /// * `amount` - Scale amount.
    pub fn new(tex1: ArcTexture<T>, tex2: ArcTexture<T>, amount: ArcTexture<Float>) -> Self {
        Self {
            tex1: tex1.clone(),
            tex2: tex2.clone(),
            amount: amount.clone(),
        }
    }
}

impl<T> Texture<T> for MixTexture<T>
where
    T: Copy + Clone + Add<Output = T>,
    Float: Mul<T, Output = T>,
{
    /// Evaluate the texture at surface interaction.
    ///
    /// * `si` - Surface interaction.
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        let t1 = self.tex1.evaluate(si);
        let t2 = self.tex2.evaluate(si);
        let amt = self.amount.evaluate(si);
        (1.0 - amt) * t1 + amt * t2
    }
}
