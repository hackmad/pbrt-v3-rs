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
    /// * `si` - The surface interaction.
    fn map(&self, si: &SurfaceInteraction) -> TextureMap2DResult {
        // Compute texture differentials for 2D identity mapping.
        let dstdx = Vector2f::new(self.su * si.dudx, self.sv * si.dvdx);
        let dstdy = Vector2f::new(self.su * si.dudy, self.sv * si.dvdy);
        let p = Point2f::new(self.su * si.uv[0] + self.du, self.sv * si.uv[1] + self.dv);
        TextureMap2DResult::new(p, dstdx, dstdy)
    }
}
