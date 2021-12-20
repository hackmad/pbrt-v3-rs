//! Scale Texture

use super::*;
use core::geometry::*;
use core::pbrt::*;
use core::spectrum::*;
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
            tex1: Arc::clone(&tex1),
            tex2: Arc::clone(&tex2),
        }
    }
}

impl<T> Texture<T> for ScaleTexture<T>
where
    T: Copy + Mul<Output = T>,
{
    /// Evaluate the texture at surface interaction.
    ///
    /// * `hit` - Surface interaction hit.
    /// * `uv`  - Surface interaction uv.
    /// * `der` - Surface interaction derivatives.
    fn evaluate(&self, hit: &Hit, uv: &Point2f, der: &Derivatives) -> T {
        self.tex1.evaluate(hit, uv, der) * self.tex2.evaluate(hit, uv, der)
    }
}

macro_rules! from_params {
    ($t: ty, $get_texture_or_else_func: ident) => {
        impl From<(&TextureParams, ArcTransform)> for ScaleTexture<$t> {
            /// Create a `ScaleTexture<$t>` from given parameter set and
            /// transformation from texture space to world space.
            ///
            /// * `p` - Tuple containing texture parameters and texture space
            ///         to world space transform.
            fn from(p: (&TextureParams, ArcTransform)) -> Self {
                let (tp, _tex2world) = p;

                let tex1 = tp.$get_texture_or_else_func("tex1", 1.0.into(), |v| {
                    Arc::new(ConstantTexture::new(v))
                });

                let tex2 = tp.$get_texture_or_else_func("tex2", 1.0.into(), |v| {
                    Arc::new(ConstantTexture::new(v))
                });

                Self::new(tex1, tex2)
            }
        }
    };
}
from_params!(Float, get_float_texture_or_else);
from_params!(Spectrum, get_spectrum_texture_or_else);
