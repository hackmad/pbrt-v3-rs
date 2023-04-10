//! 3D Checkerboard

use super::*;
use core::geometry::*;
use core::interaction::*;
use core::pbrt::*;
use core::spectrum::*;

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
        Self { tex1, tex2, mapping }
    }
}

impl<T> Texture<T> for CheckerboardTexture3D<T>
where
    T: Copy,
{
    /// Evaluate the texture at surface interaction.
    ///
    /// * `hit` - Surface interaction hit.
    /// * `uv`  - Surface interaction uv.
    /// * `der` - Surface interaction derivatives.
    fn evaluate(&self, hit: &Hit, uv: &Point2f, der: &Derivatives) -> T {
        // Get the (s, t) mapping for the intersection.
        let TextureMap3DResult { p, .. } = self.mapping.map(hit, uv, der);
        if (p.x.floor() as Int + p.y.floor() as Int + p.z.floor() as Int) % 2 == 0 {
            self.tex1.evaluate(hit, uv, der)
        } else {
            self.tex2.evaluate(hit, uv, der)
        }
    }
}

macro_rules! from_params {
    ($t: ty, $get_texture_or_else_func: ident) => {
        impl From<(&TextureParams, ArcTransform)> for CheckerboardTexture3D<$t> {
            /// Create a `CheckerboardTexture3D<$t>` from given parameter set and transformation from texture space to
            /// world space.
            ///
            /// * `p` - Tuple containing texture parameters and texture space to world space transform.
            fn from(p: (&TextureParams, ArcTransform)) -> Self {
                let (tp, tex2world) = p;

                // Check texture dimensions.
                let dim = tp.find_int("dimension", 3);
                if dim != 3 {
                    panic!("Cannot create CheckerboardTexture3D for dim = {}", dim);
                }
                // Get textures.
                let tex1 = tp.$get_texture_or_else_func("tex1", 1.0.into(), |v| Arc::new(ConstantTexture::new(v)));

                let tex2 = tp.$get_texture_or_else_func("tex2", 0.0.into(), |v| Arc::new(ConstantTexture::new(v)));

                // Initialize 3D texture mapping `map` from `tex2world`.
                let map = Arc::new(IdentityMapping3D::new(tex2world));

                Self::new(tex1, tex2, map)
            }
        }
    };
}
from_params!(Float, get_float_texture_or_else);
from_params!(Spectrum, get_spectrum_texture_or_else);
