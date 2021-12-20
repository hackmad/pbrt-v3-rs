//! 3D Identity Mapping

use super::*;

/// Implements 3D identity mapping by simply transforming the hit point and
/// partials from world space to texture space.
pub struct IdentityMapping3D {
    /// Transformation from world space to texture space.
    world_to_texture: ArcTransform,
}

impl IdentityMapping3D {
    /// Create a new `IdentityMapping3D`.
    ///
    /// * `world_to_texture` - Transformation from world space to texture space.
    pub fn new(world_to_texture: ArcTransform) -> Self {
        Self { world_to_texture }
    }
}

impl TextureMapping3D for IdentityMapping3D {
    /// Returns the (s, t) texture coordinates and partial derivitives.
    ///
    /// * `hit` - Surface interaction hit.
    /// * `uv`  - Surface interaction uv.
    /// * `der` - Surface interaction derivatives.
    fn map(&self, hit: &Hit, _uv: &Point2f, der: &Derivatives) -> TextureMap3DResult {
        let dpdx = self.world_to_texture.transform_vector(&der.dpdx);
        let dpdy = self.world_to_texture.transform_vector(&der.dpdy);
        let p = self.world_to_texture.transform_point(&hit.p);
        TextureMap3DResult::new(p, dpdx, dpdy)
    }
}
