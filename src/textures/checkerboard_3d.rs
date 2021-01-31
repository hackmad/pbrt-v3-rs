//! 3D Checkerboard

#![allow(dead_code)]
use super::TextureProps;
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::spectrum::*;
use crate::core::texture::*;
use crate::textures::*;

/// Implements a checkerboard texture via a 3D mapping.
#[derive(Clone)]
pub struct CheckerboardTexture3D<T> {
    /// First texture.
    tex1: ArcTexture<T>,

    /// Second texture.
    tex2: ArcTexture<T>,

    /// 3D mapping.
    mapping: ArcTextureMapping3D,
}

impl<T> CheckerboardTexture3D<T> {
    /// Create a new `CheckerboardTexture3D<T>`.
    ///
    /// * `tex1`      - The first texture.
    /// * `tex2`      - The second texture.
    /// * `mapping`   - The 3D mapping.
    pub fn new(tex1: ArcTexture<T>, tex2: ArcTexture<T>, mapping: ArcTextureMapping3D) -> Self {
        Self {
            tex1: tex1.clone(),
            tex2: tex2.clone(),
            mapping: mapping.clone(),
        }
    }
}

impl<T> Texture<T> for CheckerboardTexture3D<T>
where
    T: Copy,
{
    /// Evaluate the texture at surface interaction.
    ///
    /// * `si` - Surface interaction.
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        // Get the (s, t) mapping for the intersection.
        let TextureMap3DResult { p, .. } = self.mapping.map(si);
        if (p.x.floor() as Int + p.y.floor() as Int + p.z.floor() as Int) % 2 == 0 {
            self.tex1.evaluate(si)
        } else {
            self.tex2.evaluate(si)
        }
    }
}

macro_rules! from_params {
    ($t: ty, $get_texture_or_else_func: ident) => {
        impl From<&mut TextureProps> for CheckerboardTexture3D<$t> {
            /// Create a `CheckerboardTexture3D<$t>` from given parameter set and
            /// transformation from texture space to world space.
            ///
            /// * `props` - Texture creation properties.
            fn from(props: &mut TextureProps) -> Self {
                // Check texture dimensions.
                let dim = props.tp.find_int("dimension", 3);
                if dim != 3 {
                    panic!("Cannot create CheckerboardTexture3D for dim = {}", dim);
                }
                // Get textures.
                let tex1 = props
                    .tp
                    .$get_texture_or_else_func("tex1", Arc::new(ConstantTexture::new(1.0.into())));
                let tex2 = props
                    .tp
                    .$get_texture_or_else_func("tex2", Arc::new(ConstantTexture::new(0.0.into())));
                // Initialize 3D texture mapping `map` from `tex2world`.
                let map = Arc::new(IdentityMapping3D::new(props.tex2world));
                Self::new(tex1, tex2, map)
            }
        }
    };
}
from_params!(Float, get_float_texture_or_else);
from_params!(Spectrum, get_spectrum_texture_or_else);
