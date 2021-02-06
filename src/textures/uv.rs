//! UV Texture

#![allow(dead_code)]
use super::get_texture_mapping;
use crate::core::geometry::*;
use crate::core::paramset::*;
use crate::core::spectrum::*;
use crate::core::texture::*;

/// Implements a texture that converts the surface's (u, v) coordinates
/// into red and green components of a `Spectrum`.
#[derive(Clone)]
pub struct UVTexture {
    /// 2D mapping.
    mapping: ArcTextureMapping2D,
}

impl UVTexture {
    /// Create a new `UVTexture<T>`.
    ///
    /// * `mapping` - The 2D mapping.
    pub fn new(mapping: ArcTextureMapping2D) -> Self {
        Self {
            mapping: mapping.clone(),
        }
    }
}

impl Texture<Spectrum> for UVTexture {
    /// Evaluate the texture at surface interaction.
    ///
    /// * `si` - Surface interaction.
    fn evaluate(&self, si: &SurfaceInteraction) -> Spectrum {
        // Get the (s, t) mapping for the intersection.
        let TextureMap2DResult { p: st, .. } = self.mapping.map(si);

        let rgb = [st[0] - st[0].floor(), st[1] - st[1].floor(), 0.0];
        Spectrum::from_rgb(&rgb, None)
    }
}

impl From<(&TextureParams, &Transform)> for UVTexture {
    /// Create a `UVTexture` from given parameter set and
    /// transformation from texture space to world space.
    ///
    /// * `p` - Tuple containing texture parameters and texture space
    ///         to world space transform.
    fn from(p: (&TextureParams, &Transform)) -> Self {
        let (tp, tex2world) = p;

        // Initialize 2D texture mapping `map` from `tp`.
        let map = get_texture_mapping(tp, tex2world);
        Self::new(map)
    }
}
