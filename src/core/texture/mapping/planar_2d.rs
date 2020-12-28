//! 2D Planar Mapping

#![allow(dead_code)]
use super::*;
use crate::core::pbrt::*;

/// Implements 2D planar mapping by projecting points onto a plane.
pub struct PlanarMapping2D {
    /// Direction in `s` direction.
    vs: Vector3f,

    /// Direction in `t` direction (not parallel to `vs`).
    vt: Vector3f,

    /// Offset in `s` direction.
    ds: Float,

    /// Offset in `t` direction.
    dt: Float,
}

impl PlanarMapping2D {
    /// Create a new `PlanarMapping2D`.
    ///
    /// * `vs` - Direction in `s` direction.
    /// * `vt` - Direction in `t` direction (not parallel to `vs`).
    /// * `ds` - Offset in `s` direction.
    /// * `dt` - Offset in `t` direction.
    pub fn new(vs: Vector3f, vt: Vector3f, ds: Float, dt: Float) -> Self {
        Self { vs, vt, ds, dt }
    }
}

impl TextureMapping2D for PlanarMapping2D {
    /// Returns the (s, t) texture coordinates and texture differentials.
    ///
    /// * `si` - The surface interaction.
    fn map(&self, si: &SurfaceInteraction) -> TextureMap2DResult {
        let vec = Vector3f::from(si.hit.p);
        let dstdx = Vector2f(si.dpdx.dot(&self.vs), si.dpdx.dot(&self.vt));
        let dstdy = Vector2f(si.dpdy.dot(&self.vs), si.dpdy.dot(&self.vt));
        let p = Point2f::from(self.ds + vec.dot(&self.vs), dt + vec.dot(&self.vt));

        TextureMap2DResult::new(p, dstdx, dstdy)
    }
}
