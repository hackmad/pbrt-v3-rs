//! Windy Waves Texture

#[allow(dead_code)]
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::texture::*;
use std::marker::PhantomData;

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
            mapping: mapping.clone(),
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
    /// * `si` - Surface interaction.
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        // Get the (s, t) mapping for the intersection.
        let TextureMap3DResult { p, dpdx, dpdy } = self.mapping.map(si);
        let wind_strength = fbm(&(0.1 * p), &(0.1 * dpdx), &(0.1 * dpdy), 0.5, 3);
        let wave_height = fbm(&p, &dpdx, &dpdy, 0.5, 6);
        (abs(wind_strength) * wave_height).into()
    }
}
