//! FBm Texture

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::texture::*;
use std::marker::PhantomData;

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
            mapping: mapping.clone(),
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
    /// * `si` - Surface interaction.
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        // Get the (s, t) mapping for the intersection.
        let TextureMap3DResult { p, dpdx, dpdy } = self.mapping.map(si);
        fbm(&p, &dpdx, &dpdy, self.omega, self.octaves).into()
    }
}
