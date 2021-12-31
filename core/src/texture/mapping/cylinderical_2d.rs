//! 2D Cylinderical Mapping

use super::*;
use crate::pbrt::*;

/// Implements 2D cylinderical mapping.
pub struct CylindericalMapping2D {
    /// Transformation from world space to texture space.
    world_to_texture: Transform,
}

impl CylindericalMapping2D {
    /// Create a new `CylindericalMapping2D`.
    ///
    /// * `world_to_texture` - Transformation from world space to texture space.
    pub fn new(world_to_texture: Transform) -> Self {
        Self { world_to_texture }
    }

    /// Project a single point away from central axis of cylinder up to the surface.
    ///
    /// * `p` - The point.
    fn cylinder(&self, p: &Point3f) -> Point2f {
        let vec =
            (self.world_to_texture.transform_point(p) - Point3f::new(0.0, 0.0, 0.0)).normalize();
        Point2f::new((PI + atan2(vec.y, vec.x)) * INV_TWO_PI, vec.z)
    }
}

impl TextureMapping2D for CylindericalMapping2D {
    /// Returns the (s, t) texture coordinates and texture differentials.
    ///
    /// * `hit` - Surface interaction hit.
    /// * `uv`  - Surface interaction uv.
    /// * `der` - Surface interaction derivatives.
    fn map(&self, hit: &Hit, _uv: &Point2f, der: &Derivatives) -> TextureMap2DResult {
        let st = self.cylinder(&hit.p);

        // Compute texture coordinate differentials for cylinder (u, v) mapping.
        let delta = 0.01;

        let st_delta_x = self.cylinder(&(hit.p + delta * der.dpdx));
        let mut dstdx = (st_delta_x - st) / delta;
        if dstdx[1] > 0.5 {
            dstdx[1] = 1.0 - dstdx[1];
        } else if dstdx[1] < -0.5 {
            dstdx[1] = -(dstdx[1] + 1.0);
        }

        let st_delta_y = self.cylinder(&(hit.p + delta * der.dpdy));
        let mut dstdy = (st_delta_y - st) / delta;
        if dstdy[1] > 0.5 {
            dstdy[1] = 1.0 - dstdy[1];
        } else if dstdy[1] < -0.5 {
            dstdy[1] = -(dstdy[1] + 1.0);
        }

        TextureMap2DResult::new(st, dstdx, dstdy)
    }
}
