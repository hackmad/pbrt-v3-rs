//! 2D Checkerboard

#[allow(dead_code)]
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::texture::*;
use std::ops::{Add, Mul};

#[derive(Clone, PartialEq)]
#[repr(C)]
pub enum AAMethod {
    None = 0,
    ClosedForm = 1,
}

/// Implements a checkerboard texture via a 2D mapping.
#[derive(Clone)]
pub struct CheckerboardTexture2D<T> {
    /// First texture.
    tex1: ArcTexture<T>,

    /// Second texture.
    tex2: ArcTexture<T>,

    /// 2D mapping.
    mapping: ArcTextureMapping2D,

    /// Antialiasing method.
    aa_method: AAMethod,
}

impl<T> CheckerboardTexture2D<T> {
    /// Create a new `CheckerboardTexture2D<T>`.
    ///
    /// * `tex1`      - The first texture.
    /// * `tex2`      - The second texture.
    /// * `mapping`   - The 2D mapping.
    /// * `aa_method` - The antialiasing method.
    pub fn new(
        tex1: ArcTexture<T>,
        tex2: ArcTexture<T>,
        mapping: ArcTextureMapping2D,
        aa_method: AAMethod,
    ) -> Self {
        Self {
            tex1: tex1.clone(),
            tex2: tex2.clone(),
            mapping: mapping.clone(),
            aa_method,
        }
    }
}

impl<T> Texture<T> for CheckerboardTexture2D<T>
where
    T: Copy + Clone + Add<Output = T> + Mul<Float, Output = T>,
{
    /// Evaluate the texture at surface interaction.
    ///
    /// * `si` - Surface interaction.
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        // Get the (s, t) mapping for the intersection.
        let TextureMap2DResult {
            p: st,
            dstdx,
            dstdy,
        } = self.mapping.map(si);

        if self.aa_method == AAMethod::None {
            // Point sample `Checkerboard2DTexture2D`.
            if (st[0].floor() as Int + st[1].floor() as Int) % 2 == 0 {
                return self.tex1.evaluate(si);
            }
            return self.tex2.evaluate(si);
        } else {
            // Compute closed-form box-filtered `Checkerboard2DTexture2D` value.

            // Evaluate single check if filter is entirely inside one of them.
            let ds = max(abs(dstdx[0]), abs(dstdy[0]));
            let dt = max(abs(dstdx[1]), abs(dstdy[1]));

            let s0 = st[0] - ds;
            let s1 = st[0] + ds;
            let t0 = st[1] - dt;
            let t1 = st[1] + dt;
            if s0.floor() == s1.floor() && t0.floor() == t1.floor() {
                // Point sample `Checkerboard2DTexture2D`.
                if (st[0].floor() as Int + st[1].floor() as Int) % 2 == 0 {
                    return self.tex1.evaluate(si);
                }
                return self.tex2.evaluate(si);
            }

            // Apply box filter to checkerboard region.
            let bump_int = |x: Float| -> Int {
                (x / 2.0).floor() as Int + 2 * max(x / 2.0 - (x / 2.0).floor() - 0.5, 0.0) as Int
            };

            let sint = (bump_int(s1) - bump_int(s0)) as Float / (2.0 * ds);
            let tint = (bump_int(t1) - bump_int(t0)) as Float / (2.0 * dt);
            let area2 = if ds > 1.0 || dt > 1.0 {
                0.5
            } else {
                sint + tint - 2.0 * sint * tint
            };

            self.tex1.evaluate(si) * (1.0 - area2) + self.tex2.evaluate(si) * area2
        }
    }
}
