//! 2D Spherical Mapping

use super::*;
use crate::pbrt::*;

/// Implements 2D spherical mapping.
pub struct SphericalMapping2D {
    /// Transformation from world space to texture space.
    world_to_texture: Transform,
}

impl SphericalMapping2D {
    /// Create a new `SphericalMapping2D`.
    ///
    /// * `world_to_texture` - Transformation from world space to texture space.
    pub fn new(world_to_texture: Transform) -> Self {
        Self { world_to_texture }
    }

    /// Project a single point along the vector from the sphere's center through
    /// it up to the surface.
    ///
    /// * `p` - The point.
    fn sphere(&self, p: &Point3f) -> Point2f {
        let vec =
            (self.world_to_texture.transform_point(p) - Point3f::new(0.0, 0.0, 0.0)).normalize();
        let theta = spherical_theta(&vec);
        let phi = spherical_phi(&vec);
        Point2f::new(theta * INV_PI, phi * INV_TWO_PI)
    }
}

impl TextureMapping2D for SphericalMapping2D {
    /// Returns the (s, t) texture coordinates and texture differentials.
    ///
    /// * `si` - The surface interaction.
    fn map(&self, si: &SurfaceInteraction) -> TextureMap2DResult {
        let st = self.sphere(&si.hit.p);

        // Compute texture coordinate differentials for sphere (u, v) mapping.
        let delta = 0.1;

        let st_delta_x = self.sphere(&(si.hit.p + delta * si.dpdx));
        let mut dstdx = (st_delta_x - st) / delta;

        let st_delta_y = self.sphere(&(si.hit.p + delta * si.dpdy));
        let mut dstdy = (st_delta_y - st) / delta;

        // Handle sphere mapping discontinuity for coordinate differentials.
        if dstdx[1] > 0.5 {
            dstdx[1] = 1.0 - dstdx[1];
        } else if dstdx[1] < -0.5 {
            dstdx[1] = -(dstdx[1] + 1.0);
        }

        if dstdy[1] > 0.5 {
            dstdy[1] = 1.0 - dstdy[1];
        } else if dstdy[1] < -0.5 {
            dstdy[1] = -(dstdy[1] + 1.0);
        }

        TextureMap2DResult::new(st, dstdx, dstdy)
    }
}
