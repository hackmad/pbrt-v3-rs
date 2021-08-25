//! Texture Mapping

use crate::geometry::*;
use std::sync::Arc;

/// Stores the result of 2D texture mapping.
#[derive(Copy, Clone, Default)]
pub struct TextureMap2DResult {
    /// Texture mapped point.
    pub p: Point2f,

    /// Texture differential of (s, t) with respect to x, (∂s/∂x, ∂t/∂x).
    pub dstdx: Vector2f,

    /// Texture differential of (s, t) with respect to y, (∂s/∂y, ∂t/∂y).
    pub dstdy: Vector2f,
}

impl TextureMap2DResult {
    /// Create a new `TextureMap2DResult`.
    ///
    /// * `p`      - Texture mapped point.
    /// * `dstdx`  - Texture differential of (s, t) with respect to x, (∂s/∂x, ∂t/∂x).
    /// * `dxtdy`  - Texture differential of (s, t) with respect to y, (∂s/∂y, ∂t/∂y).
    pub fn new(p: Point2f, dstdx: Vector2f, dstdy: Vector2f) -> Self {
        Self { p, dstdx, dstdy }
    }
}

/// Interface for 2D texture mapping.
pub trait TextureMapping2D {
    /// Returns the (s, t) texture coordinates and texture differentials.
    ///
    /// * `si` - The surface interaction.
    fn map(&self, si: &SurfaceInteraction) -> TextureMap2DResult;
}

/// Atomic reference counted `TextureMapping2D`.
pub type ArcTextureMapping2D = Arc<dyn TextureMapping2D + Send + Sync>;

/// Stores the result of 3D texture mapping.
#[derive(Copy, Clone, Default)]
pub struct TextureMap3DResult {
    /// Texture mapped point.
    pub p: Point3f,

    /// Partial derivatives of position p to x, ∂p/∂x.
    pub dpdx: Vector3f,

    /// Partial derivatives of position p to y, ∂p/∂y.
    pub dpdy: Vector3f,
}

impl TextureMap3DResult {
    /// Create a new `TextureMap2DResult`.
    ///
    /// * `p`     - Texture mapped point.
    /// * `dpdx`  - Partial derivatives of position p to x, ∂p/∂x.
    /// * `dpdy`  - Partial derivatives of position p to y, ∂p/∂y.
    pub fn new(p: Point3f, dpdx: Vector3f, dpdy: Vector3f) -> Self {
        Self { p, dpdx, dpdy }
    }
}
/// Interface for 3D texture mapping.
pub trait TextureMapping3D {
    /// Returns the (s, t) texture coordinates and partial derivitives.
    ///
    /// * `si` - The surface interaction.
    fn map(&self, si: &SurfaceInteraction) -> TextureMap3DResult;
}

/// Atomic reference counted `TextureMapping3D`.
pub type ArcTextureMapping3D = Arc<dyn TextureMapping3D + Send + Sync>;

mod cylinderical_2d;
mod identity_3d;
mod planar_2d;
mod spherical_2d;
mod uv_2d;

// Re-export
pub use cylinderical_2d::*;
pub use identity_3d::*;
pub use planar_2d::*;
pub use spherical_2d::*;
pub use uv_2d::*;
