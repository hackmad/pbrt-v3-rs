//! Polka Dots

use super::get_texture_mapping;
use super::*;
use core::geometry::*;
use core::pbrt::*;
use core::spectrum::*;
use std::sync::Arc;

/// Implements a random polka dot texture via a 2D mapping.
#[derive(Clone)]
pub struct DotsTexture<T> {
    /// Texture for exterior of dot pattern.
    outside_dot: ArcTexture<T>,

    /// Texture for dot pattern.
    inside_dot: ArcTexture<T>,

    /// 2D mapping.
    mapping: ArcTextureMapping2D,
}

impl<T> DotsTexture<T> {
    /// Create a new `DotsTexture<T>`.
    ///
    /// * `outside_dot` - Texture for exterior of dot pattern.
    /// * `inside_dot`  - Texture for dot pattern.
    /// * `mapping`     - The 2D mapping.
    pub fn new(
        outside_dot: ArcTexture<T>,
        inside_dot: ArcTexture<T>,
        mapping: ArcTextureMapping2D,
    ) -> Self {
        Self {
            outside_dot: Arc::clone(&outside_dot),
            inside_dot: Arc::clone(&inside_dot),
            mapping: Arc::clone(&mapping),
        }
    }
}

impl<T> Texture<T> for DotsTexture<T>
where
    T: Copy,
{
    /// Evaluate the texture at surface interaction.
    ///
    /// * `si` - Surface interaction.
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        // Get the (s, t) mapping for the intersection.
        let TextureMap2DResult { p: st, .. } = self.mapping.map(si);

        let s_cell = (st[0] + 0.5).floor();
        let t_cell = (st[1] + 0.5).floor();

        // Return _insideDot_ result if point is inside dot
        if noise(Point3f::new(s_cell + 0.5, t_cell + 0.5, 0.0)) > 0.0 {
            let radius = 0.35;
            let max_shift = 0.5 - radius;

            let s_center =
                s_cell + max_shift * noise(Point3f::new(s_cell + 1.5, t_cell + 2.8, 0.0));
            let t_center =
                t_cell + max_shift * noise(Point3f::new(s_cell + 4.5, t_cell + 9.8, 0.0));

            let dst = st - Point2f::new(s_center, t_center);
            if dst.length_squared() < radius * radius {
                return self.inside_dot.evaluate(si);
            }
        }
        self.outside_dot.evaluate(si)
    }
}

macro_rules! from_params {
    ($t: ty, $get_texture_or_else_func: ident) => {
        impl From<(&TextureParams, ArcTransform)> for DotsTexture<$t> {
            /// Create a `DotsTexture<$t>` from given parameter set and
            /// transformation from texture space to world space.
            ///
            /// * `p` - Tuple containing texture parameters and texture space
            ///         to world space transform.
            fn from(p: (&TextureParams, ArcTransform)) -> Self {
                let (tp, tex2world) = p;

                // Initialize 2D texture mapping `map` from `tp`.
                let map = get_texture_mapping(tp, tex2world);

                let inside = tp.$get_texture_or_else_func("inside", 1.0.into(), |v| {
                    Arc::new(ConstantTexture::new(v))
                });

                let outside = tp.$get_texture_or_else_func("outside", 0.0.into(), |v| {
                    Arc::new(ConstantTexture::new(v))
                });

                Self::new(inside, outside, map)
            }
        }
    };
}
from_params!(Float, get_float_texture_or_else);
from_params!(Spectrum, get_spectrum_texture_or_else);
