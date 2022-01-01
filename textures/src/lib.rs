//! Textures

use core::geometry::{ArcTransform, Vector3f};
use core::paramset::TextureParams;
use core::texture::*;
use std::sync::Arc;

#[macro_use]
extern crate log;

mod bilerp;
mod checkerboard_2d;
mod checkerboard_3d;
mod constant;
mod dots;
mod fbm;
mod imagemap;
mod marble;
mod mix;
mod scale;
mod uv;
mod windy;
mod wrinkled;

// Re-export
pub use bilerp::*;
pub use checkerboard_2d::*;
pub use checkerboard_3d::*;
pub use constant::*;
pub use dots::*;
pub use fbm::*;
pub use imagemap::*;
pub use marble::*;
pub use mix::*;
pub use scale::*;
pub use uv::*;
pub use windy::*;
pub use wrinkled::*;

/// Returns a 2D texture mapping reference from the texture parameters.
///
/// * `tp`        - Texture parameters.
/// * `tex2world` - Texture space to world space transform.
fn get_texture_mapping(tp: &TextureParams, tex2world: ArcTransform) -> ArcTextureMapping2D {
    let map_type = tp.find_string("mapping", String::from("uv"));
    match &map_type[..] {
        "uv" => {
            let su = tp.find_float("uscale", 1.0);
            let sv = tp.find_float("vscale", 1.0);
            let du = tp.find_float("udelta", 0.0);
            let dv = tp.find_float("vdelta", 0.0);
            Arc::new(UVMapping2D::new(su, sv, du, dv))
        }
        "spherical" => Arc::new(SphericalMapping2D::new(tex2world.inverse())),
        "cylinderical" => Arc::new(CylindericalMapping2D::new(tex2world.inverse())),
        "planar" => Arc::new(PlanarMapping2D::new(
            tp.find_vector3f("v1", Vector3f::new(1.0, 0.0, 0.0)),
            tp.find_vector3f("v2", Vector3f::new(0.0, 1.0, 0.0)),
            tp.find_float("udelta", 0.0),
            tp.find_float("vdelta", 0.0),
        )),
        mt => {
            warn!("Error 2D texture mapping '{}' unknown", mt);
            Arc::new(UVMapping2D::default())
        }
    }
}
