//! Windy Waves Texture

use core::geometry::*;
use core::interaction::*;
use core::paramset::*;
use core::pbrt::*;
use core::texture::*;
use std::marker::PhantomData;
use std::sync::Arc;

/// Implements windy waves texture via a 3D mapping.
#[derive(Clone)]
pub struct WindyTexture<T> {
    /// 2D mapping.
    mapping: ArcTextureMapping3D,

    /// Compiler hint.
    _marker: PhantomData<T>,
}

impl<T> WindyTexture<T> {
    /// Create a new `WindyTexture<T>`.
    ///
    /// * `mapping` - The 3D mapping.
    pub fn new(mapping: ArcTextureMapping3D) -> Self {
        Self {
            mapping,
            _marker: PhantomData,
        }
    }
}

impl<T> Texture<T> for WindyTexture<T>
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
        let wind_strength = fbm(&(0.1 * p), &(0.1 * dpdx), &(0.1 * dpdy), 0.5, 3);
        let wave_height = fbm(&p, &dpdx, &dpdy, 0.5, 6);
        (abs(wind_strength) * wave_height).into()
    }
}

impl<T> From<(&TextureParams, ArcTransform)> for WindyTexture<T> {
    /// Create a `WindyTexture<T>` from given parameter set and transformation from texture space to world space.
    ///
    /// * `p` - Tuple containing texture parameters and texture space to world space transform.
    fn from(p: (&TextureParams, ArcTransform)) -> Self {
        let (_tp, tex2world) = p;
        let map = Arc::new(IdentityMapping3D::new(tex2world));
        Self::new(map)
    }
}
