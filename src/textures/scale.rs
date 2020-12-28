//! Scale Texture

#[allow(dead_code)]
use crate::core::geometry::*;
use crate::core::texture::*;
use std::ops::Mul;

/// Implements a texture that returns the product of 2 textures.
#[derive(Clone)]
pub struct ScaleTexture<T> {
    /// First texture.
    tex1: ArcTexture<T>,

    /// Second texture.
    tex2: ArcTexture<T>,
}

impl<T> ScaleTexture<T> {
    /// Create a new `ScaleTexture`.
    ///
    /// * `tex1` - The first texture.
    /// * `tex2` - The second texture.
    pub fn new(tex1: ArcTexture<T>, tex2: ArcTexture<T>) -> Self {
        Self {
            tex1: tex1.clone(),
            tex2: tex2.clone(),
        }
    }
}

impl<T> Texture<T> for ScaleTexture<T>
where
    T: Copy + Mul<Output = T>,
{
    /// Evaluate the texture at surface interaction.
    ///
    /// * `si` - Surface interaction.
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        self.tex1.evaluate(si) * self.tex2.evaluate(si)
    }
}
