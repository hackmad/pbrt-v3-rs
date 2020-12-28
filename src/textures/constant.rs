//! Constant Texture

#[allow(dead_code)]
use crate::core::geometry::*;
use crate::core::texture::*;

/// Implements a texture that returns the same value everywhere.
pub struct ConstantTexture<T> {
    /// The texture value.
    value: T,
}

impl<T: Copy> ConstantTexture<T> {
    /// Create a new `ConstantTexture<T>`.
    ///
    /// * `value` - The texture value.
    pub fn new(value: T) -> Self {
        Self { value }
    }
}

impl<T: Copy> Texture<T> for ConstantTexture<T> {
    /// Evaluate the texture at surface interaction.
    ///
    /// * `si` - Surface interaction.
    fn evaluate(&self, _si: &SurfaceInteraction) -> T {
        self.value
    }
}
