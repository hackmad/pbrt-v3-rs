//! Marble Texture

use core::geometry::*;
use core::interaction::*;
use core::paramset::*;
use core::pbrt::*;
use core::spectrum::*;
use core::texture::*;
use std::sync::Arc;

/// Implements marble texture via a 3D mapping.
#[derive(Clone)]
pub struct MarbleTexture {
    /// 3D mapping.
    mapping: ArcTextureMapping3D,

    /// Smoothness falloff value in [0, 1].
    omega: Float,

    /// Maximum number of octaves of noise to use for the sum.
    octaves: usize,

    /// Scale factor.
    scale: Float,

    /// Modulates the magnitude of perturbation.
    variation: Float,
}

impl MarbleTexture {
    /// Create a new `MarbleTexture<T>`.
    ///
    /// * `mapping`   - The 3D mapping.
    /// * `omega`     - Smoothness falloff value in [0, 1].
    /// * `octaves`   - Maximum number of octaves of noise to use for the sum.
    /// * `scale`     - Scale factor.
    /// * `variation` - Modulates the magnitude of perturbation.
    pub fn new(mapping: ArcTextureMapping3D, omega: Float, octaves: usize, scale: Float, variation: Float) -> Self {
        Self {
            mapping,
            omega,
            octaves,
            scale,
            variation,
        }
    }
}

impl Texture<Spectrum> for MarbleTexture {
    /// Evaluate the texture at surface interaction.
    ///
    /// * `hit` - Surface interaction hit.
    /// * `uv`  - Surface interaction uv.
    /// * `der` - Surface interaction derivatives.
    fn evaluate(&self, hit: &Hit, uv: &Point2f, der: &Derivatives) -> Spectrum {
        // Get the (s, t) mapping for the intersection.
        let TextureMap3DResult { p, dpdx, dpdy } = self.mapping.map(hit, uv, der);
        let p = p * self.scale;

        let marble =
            p.y + self.variation * fbm(&p, &(self.scale * dpdx), &(self.scale * dpdy), self.omega, self.octaves);
        let t = 0.5 + 0.5 * sin(marble);

        // Evaluate marble spline at `t`.
        let first = min(1, (t * NSEGS as Float).floor() as usize);
        let t = t * NSEGS as Float - first as Float;

        let c0 = Spectrum::from_rgb(&C[first], None);
        let c1 = Spectrum::from_rgb(&C[first + 1], None);
        let c2 = Spectrum::from_rgb(&C[first + 2], None);
        let c3 = Spectrum::from_rgb(&C[first + 3], None);

        // Bezier spline evaluated with de Castilejau's algorithm.
        let s0 = (1.0 - t) * c0 + t * c1;
        let s1 = (1.0 - t) * c1 + t * c2;
        let s2 = (1.0 - t) * c2 + t * c3;

        let s0 = (1.0 - t) * s0 + t * s1;
        let s1 = (1.0 - t) * s1 + t * s2;

        // Extra scale of 1.5 to increase variation among colors.
        1.5 * ((1.0 - t) * s0 + t * s1)
    }
}

impl From<(&TextureParams, ArcTransform)> for MarbleTexture {
    /// Create a `MarbleTexture<T>` from given parameter set and transformation from texture space to world space.
    ///
    /// * `p` - Tuple containing texture parameters and texture space to world space transform.
    fn from(p: (&TextureParams, ArcTransform)) -> Self {
        let (tp, tex2world) = p;
        let map = Arc::new(IdentityMapping3D::new(tex2world));
        Self::new(
            map,
            tp.find_float("roughness", 0.5),
            tp.find_int("octaves", 8_i32) as usize,
            tp.find_float("scale", 1.0),
            tp.find_float("variation", 0.2),
        )
    }
}

/// The marble spline coefficients.
const C: [[Float; 3]; 9] = [
    [0.58, 0.58, 0.6],
    [0.58, 0.58, 0.6],
    [0.58, 0.58, 0.6],
    [0.5, 0.5, 0.5],
    [0.6, 0.59, 0.58],
    [0.58, 0.58, 0.6],
    [0.58, 0.58, 0.6],
    [0.2, 0.2, 0.33],
    [0.58, 0.58, 0.6],
];

/// Number of spline segments.
const NSEGS: usize = 6;
