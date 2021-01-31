//! Textures

use crate::core::geometry::{Transform, Vector3f};
use crate::core::paramset::TextureParams;
use crate::core::texture::*;
use std::sync::Arc;

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

/// Stores properties for texture creation.
#[derive(Clone)]
pub struct TextureProps {
    /// Texture parameter set.
    pub tp: TextureParams,

    /// Transformation from texture space to world space.
    pub tex2world: Transform,
}

/// Returns a 2D texture mapping reference from the texture parameters.
///
/// * `props` - Texture creation properties.
fn get_texture_mapping(props: &mut TextureProps) -> ArcTextureMapping2D {
    let map_type = props.tp.find_string("mapping", String::from("uv"));
    match &map_type[..] {
        "uv" => {
            let su = props.tp.find_float("uscale", 1.0);
            let sv = props.tp.find_float("vscale", 1.0);
            let du = props.tp.find_float("udelta", 0.0);
            let dv = props.tp.find_float("vdelta", 0.0);
            Arc::new(UVMapping2D::new(su, sv, du, dv))
        }
        "spherical" => Arc::new(SphericalMapping2D::new(props.tex2world.inverse())),
        "cylinderical" => Arc::new(CylindericalMapping2D::new(props.tex2world.inverse())),
        "planar" => Arc::new(PlanarMapping2D::new(
            props.tp.find_vector3f("v1", Vector3f::new(1.0, 0.0, 0.0)),
            props.tp.find_vector3f("v2", Vector3f::new(0.0, 1.0, 0.0)),
            props.tp.find_float("udelta", 0.0),
            props.tp.find_float("vdelta", 0.0),
        )),
        mt => {
            eprintln!("Error 2D texture mapping '{}' unknown", mt);
            Arc::new(UVMapping2D::default())
        }
    }
}
