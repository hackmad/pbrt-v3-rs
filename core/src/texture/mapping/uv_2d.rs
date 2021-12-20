//! 2D (u, v) Mapping

use super::*;
use crate::pbrt::*;

/// Implements 2D (u, v) mapping.
pub struct UVMapping2D {
    /// Scale `u`.
    su: Float,

    /// Scale `v`.
    sv: Float,

    /// Shift `u`.
    du: Float,

    /// Shift `v`.
    dv: Float,
}

impl UVMapping2D {
    /// Create a new `UVMapping2D` with scale and shift values.
    ///
    /// * `su` - Scale `u`.
    /// * `sv` - Scale `v`.
    /// * `du` - Shift `u`.
    /// * `dv` - Shift `v`.
    pub fn new(su: Float, sv: Float, du: Float, dv: Float) -> Self {
        Self { su, sv, du, dv }
    }
}

impl Default for UVMapping2D {
    /// Returns a default value for `UVMapping2D` with no scaling or shifting.
    fn default() -> Self {
        Self::new(1.0, 1.0, 0.0, 0.0)
    }
}

impl TextureMapping2D for UVMapping2D {
    /// Returns the (s, t) texture coordinates and texture differentials.
    ///
    /// * `hit` - Surface interaction hit.
    /// * `uv`  - Surface interaction uv.
    /// * `der` - Surface interaction derivatives.
    fn map(&self, _hit: &Hit, uv: &Point2f, der: &Derivatives) -> TextureMap2DResult {
        // Compute texture differentials for 2D identity mapping.
        let dstdx = Vector2f::new(self.su * der.dudx, self.sv * der.dvdx);
        let dstdy = Vector2f::new(self.su * der.dudy, self.sv * der.dvdy);
        let p = Point2f::new(self.su * uv[0] + self.du, self.sv * uv[1] + self.dv);
        TextureMap2DResult::new(p, dstdx, dstdy)
    }
}
