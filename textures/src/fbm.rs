//! FBm Texture

use core::geometry::*;
use core::interaction::*;
use core::paramset::*;
use core::pbrt::*;
use core::texture::*;
use std::marker::PhantomData;
use std::sync::Arc;

/// Implements FBm (Fractional Brownian motion) texture via a 3D mapping.
#[derive(Clone)]
pub struct FBmTexture<T> {
    /// 3D mapping.
    mapping: ArcTextureMapping3D,

    /// Smoothness falloff value in [0, 1].
    omega: Float,

    /// Maximum number of octaves of noise to use for the sum.
    octaves: usize,

    /// Compiler hint.
    _marker: PhantomData<T>,
}

impl<T> FBmTexture<T> {
    /// Create a new `FBmTexture<T>`.
    ///
    /// * `mapping` - The 3D mapping.
    /// * `omega`   - Smoothness falloff value in [0, 1].
    /// * `octaves` - Maximum number of octaves of noise to use for the sum.
    pub fn new(mapping: ArcTextureMapping3D, omega: Float, octaves: usize) -> Self {
        Self {
            mapping,
            omega,
            octaves,
            _marker: PhantomData,
        }
    }
}

impl<T> Texture<T> for FBmTexture<T>
where
    T: Copy + From<Float>,
{
    /// Evaluate the texture at surface interaction.
    ///
    /// * `hit` - Surface interaction hit.
    /// * `uv`  - Surface interaction uv.
    /// * `der` - Surface interaction derivatives.
    fn evaluate(&self, hit: &Hit, uv: &Point2f, der: &Derivatives) -> T {
        // Get the (s, t) mapping for the intersection.
        let TextureMap3DResult { p, dpdx, dpdy } = self.mapping.map(hit, uv, der);
        fbm(&p, &dpdx, &dpdy, self.omega, self.octaves).into()
    }
}

impl<T> From<(&TextureParams, ArcTransform)> for FBmTexture<T> {
    /// Create a `FBmTexture<T>` from given parameter set and
    /// transformation from texture space to world space.
    ///
    /// * `p` - Tuple containing texture parameters and texture space
    ///         to world space transform.
    fn from(p: (&TextureParams, ArcTransform)) -> Self {
        let (tp, tex2world) = p;
        let map = Arc::new(IdentityMapping3D::new(tex2world));
        Self::new(
            map,
            tp.find_float("roughness", 0.5),
            tp.find_int("octaves", 8_i32) as usize,
        )
    }
}
