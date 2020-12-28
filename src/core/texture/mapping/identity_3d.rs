//! 3D Identity Mapping

#[allow(dead_code)]
use super::*;

/// Implements 3D identity mapping by simply transforming the hit point and
/// partials from world space to texture space.
pub struct IdentityMapping3D {
    /// Transformation from world space to texture space.
    world_to_texture: Transform,
}

impl IdentityMapping3D {
    /// Create a new `IdentityMapping3D`.
    ///
    /// * `world_to_texture` - Transformation from world space to texture space.
    pub fn new(world_to_texture: Transform) -> Self {
        Self { world_to_texture }
    }
}

impl TextureMapping3D for IdentityMapping3D {
    /// Returns the (s, t) texture coordinates and partial derivitives.
    ///
    /// * `si` - The surface interaction.
    fn map(&self, si: &SurfaceInteraction) -> TextureMap3DResult {
        let dpdx = self.world_to_texture.transform_vector(&si.dpdx);
        let dpdy = self.world_to_texture.transform_vector(&si.dpdy);
        let p = self.world_to_texture.transform_point(&si.hit.p);
        TextureMap3DResult::new(p, dpdx, dpdy)
    }
}
