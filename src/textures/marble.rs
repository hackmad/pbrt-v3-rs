//! Marble Texture

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::spectrum::*;
use crate::core::texture::*;

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
    pub fn new(
        mapping: ArcTextureMapping3D,
        omega: Float,
        octaves: usize,
        scale: Float,
        variation: Float,
    ) -> Self {
        Self {
            mapping: mapping.clone(),
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
    /// * `si` - Surface interaction.
    fn evaluate(&self, si: &SurfaceInteraction) -> Spectrum {
        // Get the (s, t) mapping for the intersection.
        let TextureMap3DResult { p, dpdx, dpdy } = self.mapping.map(si);
        let p = p * self.scale;

        let marble = p.y
            + self.variation
                * fbm(
                    &p,
                    &(self.scale * dpdx),
                    &(self.scale * dpdy),
                    self.omega,
                    self.octaves,
                );
        let t = 0.5 + 0.5 * sin(marble);

        // Evaluate marble spline at `t`.
        let first = min(1, (t * NSEGS as Float).floor() as usize);
        let t = t * NSEGS as Float - first as Float;

        let mut c0 = Spectrum::default();
        let mut c1 = Spectrum::default();
        let mut c2 = Spectrum::default();
        let mut c3 = Spectrum::default();

        c0.from_rgb(&C[first], SpectrumType::Reflectance);
        c1.from_rgb(&C[first + 1], SpectrumType::Reflectance);
        c2.from_rgb(&C[first + 2], SpectrumType::Reflectance);
        c3.from_rgb(&C[first + 3], SpectrumType::Reflectance);

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
