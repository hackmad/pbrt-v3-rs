//! Mix Texture

use super::*;
use core::geometry::*;
use core::interaction::*;
use core::pbrt::*;
use core::spectrum::*;
use std::ops::{Add, Mul};
use std::sync::Arc;

/// Implements a texture that linearly interpolates between two textures with a third texture.
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
        Self { tex1, tex2, amount }
    }
}

impl<T> Texture<T> for MixTexture<T>
where
    T: Copy + Clone + Add<Output = T>,
    Float: Mul<T, Output = T>,
{
    /// Evaluate the texture at surface interaction.
    ///
    /// * `hit` - Surface interaction hit.
    /// * `uv`  - Surface interaction uv.
    /// * `der` - Surface interaction derivatives.
    fn evaluate(&self, hit: &Hit, uv: &Point2f, der: &Derivatives) -> T {
        let t1 = self.tex1.evaluate(hit, uv, der);
        let t2 = self.tex2.evaluate(hit, uv, der);
        let amt = self.amount.evaluate(hit, uv, der);
        (1.0 - amt) * t1 + amt * t2
    }
}

macro_rules! from_params {
    ($t: ty, $get_texture_or_else_func: ident) => {
        impl From<(&TextureParams, ArcTransform)> for MixTexture<$t> {
            /// Create a `MixTexture<$t>` from given parameter set and transformation from texture space to world space.
            ///
            /// * `p` - Tuple containing texture parameters and texture space to world space transform.
            fn from(p: (&TextureParams, ArcTransform)) -> Self {
                let (tp, _tex2world) = p;

                let tex1 = tp.$get_texture_or_else_func("tex1", 0.0.into(), |v| Arc::new(ConstantTexture::new(v)));

                let tex2 = tp.$get_texture_or_else_func("tex2", 1.0.into(), |v| Arc::new(ConstantTexture::new(v)));

                let amt = tp.get_float_texture_or_else("amount", 0.5, |v| Arc::new(ConstantTexture::new(v)));

                Self::new(tex1, tex2, amt)
            }
        }
    };
}
from_params!(Float, get_float_texture_or_else);
from_params!(Spectrum, get_spectrum_texture_or_else);
