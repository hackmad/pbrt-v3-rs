//! Scale Texture

#![allow(dead_code)]
use super::TextureProps;
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::spectrum::*;
use crate::core::texture::*;
use crate::textures::*;
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

macro_rules! from_params {
    ($t: ty, $get_texture_or_else_func: ident) => {
        impl From<&mut TextureProps> for ScaleTexture<$t> {
            /// Create a `ScaleTexture<$t>` from given parameter set and
            /// transformation from texture space to world space.
            ///
            /// * `props` - Texture creation properties.
            fn from(props: &mut TextureProps) -> Self {
                let tex1 = props
                    .tp
                    .$get_texture_or_else_func("tex1", Arc::new(ConstantTexture::new(1.0.into())));
                let tex2 = props
                    .tp
                    .$get_texture_or_else_func("tex2", Arc::new(ConstantTexture::new(1.0.into())));
                Self::new(tex1, tex2)
            }
        }
    };
}
from_params!(Float, get_float_texture_or_else);
from_params!(Spectrum, get_spectrum_texture_or_else);
