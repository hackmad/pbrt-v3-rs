//! Bilinear Interpolation Texture

use super::get_texture_mapping;
use core::geometry::*;
use core::paramset::*;
use core::pbrt::*;
use core::spectrum::*;
use core::texture::*;
use std::ops::{Add, Mul};
use std::sync::Arc;

/// Implements a texture that bilinearly interpolates among four constant values.
#[derive(Clone)]
pub struct BilerpTexture<T> {
    /// Value at (0, 0).
    v00: T,

    /// Value at (0, 1).
    v01: T,

    /// Value at (1, 0).
    v10: T,

    /// Value at (1, 1).
    v11: T,

    /// 2D mapping.
    mapping: ArcTextureMapping2D,
}

impl<T: Copy> BilerpTexture<T> {
    /// Create a new `BilerpTexture<T>`.
    ///
    /// * `v00`     - Value at (0, 0).
    /// * `v01`     - Value at (0, 1).
    /// * `v10`     - Value at (1, 0).
    /// * `v11`     - Value at (1, 1).
    /// * `mapping` - 2D mapping.
    pub fn new(v00: T, v01: T, v10: T, v11: T, mapping: ArcTextureMapping2D) -> Self {
        Self {
            v00,
            v01,
            v10,
            v11,
            mapping: Arc::clone(&mapping),
        }
    }
}

impl<T> Texture<T> for BilerpTexture<T>
where
    T: Copy + Add<Output = T> + Mul<Float, Output = T>,
{
    /// Evaluate the texture at surface interaction.
    ///
    /// * `si` - Surface interaction.
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        // Get the (s, t) mapping for the intersection.
        let TextureMap2DResult { p: st, .. } = self.mapping.map(si);

        // Calculate the weights for the 4 points.
        let s00 = (1.0 - st[0]) * (1.0 - st[1]);
        let s01 = (1.0 - st[0]) * st[1];
        let s10 = st[0] * (1.0 - st[1]);
        let s11 = st[0] * st[1];

        // Return the interpolated value.
        (self.v00 * s00) + (self.v01 * s01) + (self.v10 * s10) + (self.v11 * s11)
    }
}

macro_rules! from_params {
    ($t: ty, $find_func: ident) => {
        impl From<(&TextureParams, &Transform)> for BilerpTexture<$t> {
            /// Create a `BilerpTexture<$t>` from given parameter set and
            /// transformation from texture space to world space.
            ///
            /// * `p` - Tuple containing texture parameters and texture space
            ///         to world space transform.
            fn from(p: (&TextureParams, &Transform)) -> Self {
                let (tp, tex2world) = p;

                // Initialize 2D texture mapping `map` from `tp`.
                let map = get_texture_mapping(tp, tex2world);
                Self::new(
                    tp.$find_func("v00", 0.0.into()),
                    tp.$find_func("v01", 1.0.into()),
                    tp.$find_func("v10", 0.0.into()),
                    tp.$find_func("v11", 1.0.into()),
                    map,
                )
            }
        }
    };
}
from_params!(Float, find_float);
from_params!(Spectrum, find_spectrum);
