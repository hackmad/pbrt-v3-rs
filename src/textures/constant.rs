//! Constant Texture

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::paramset::*;
use crate::core::pbrt::*;
use crate::core::spectrum::*;
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

macro_rules! from_params {
    ($t: ty, $find_func: ident) => {
        impl From<(&TextureParams, &Transform)> for ConstantTexture<$t> {
            /// Create a `ConstantTexture<$t>` from given parameter set and
            /// transformation from texture space to world space.
            ///
            /// * `p` - Tuple containing texture parameters and texture space
            ///         to world space transform.
            fn from(p: (&TextureParams, &Transform)) -> Self {
                let (tp, _tex2world) = p;
                Self::new(tp.$find_func("value", 1.0.into()))
            }
        }
    };
}
from_params!(Float, find_float);
from_params!(Spectrum, find_spectrum);
