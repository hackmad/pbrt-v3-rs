//! 2D Planar Mapping

use super::*;
use crate::pbrt::*;

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
    /// * `hit` - Surface interaction hit.
    /// * `uv`  - Surface interaction uv.
    /// * `der` - Surface interaction derivatives.
    fn map(&self, hit: &Hit, _uv: &Point2f, der: &Derivatives) -> TextureMap2DResult {
        let vec = Vector3f::from(hit.p);
        let dstdx = Vector2f::new(der.dpdx.dot(&self.vs), der.dpdx.dot(&self.vt));
        let dstdy = Vector2f::new(der.dpdy.dot(&self.vs), der.dpdy.dot(&self.vt));
        let p = Point2f::new(self.ds + vec.dot(&self.vs), self.dt + vec.dot(&self.vt));

        TextureMap2DResult::new(p, dstdx, dstdy)
    }
}
