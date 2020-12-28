//! Constant Texture

#[allow(dead_code)]
use crate::core::geometry::*;
use crate::core::texture::*;

/// Implements a texture that returns the same value everywhere.
#[derive(Clone)]
pub struct ConstantTexture<T> {
    /// The texture value.
    value: T,
}

impl<T> ConstantTexture<T> {
    /// Create a new `ConstantTexture<T>`.
    ///
    /// * `value` - The texture value.
    pub fn new(value: T) -> Self {
        Self { value }
    }
}

impl<T> Texture<T> for ConstantTexture<T>
where
    T: Copy,
{
    /// Evaluate the texture at surface interaction.
    ///
    /// * `si` - Surface interaction.
    fn evaluate(&self, _si: &SurfaceInteraction) -> T {
        self.value
    }
}
