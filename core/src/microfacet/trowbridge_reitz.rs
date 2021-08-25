//! Trowbridge-Reitz Distribution

#![allow(dead_code)]
use crate::geometry::*;
use crate::pbrt::*;
use crate::reflection::*;
use super::MicrofacetDistribution;

/// Implements the anisotropic variant of the Trowbridge-Reitz distribution.
#[derive(Copy, Clone, Default)]
pub struct TrowbridgeReitzDistribution {
    /// Indicates whether or not the visible area is sampled or not.
    sample_visible_area: bool,

    /// For microfacets oriented perpendicular to the x-axis and where
    /// α = sqrt(2) * σ and σ is the RMS slope of microfacets.
    alpha_x: Float,

    /// For microfacets oriented perpendicular to the y-axis and where
    /// α = sqrt(2) * σ and σ is the RMS slope of microfacets.
    alpha_y: Float,
}

impl TrowbridgeReitzDistribution {
    /// Create a new `TrowbridgeReitzDistribution`.
    ///
    /// * `alpha_x`             - For microfacets oriented perpendicular to the
    ///                           x-axis and where α = sqrt(2) * σ and σ is the
    ///                           RMS slope of microfacets.
    /// * `alpha_y`             - For microfacets oriented perpendicular to the
    ///                           y-axis and where α = sqrt(2) * σ and σ is the
    ///                           RMS slope of microfacets.
    /// * `sample_visible_area` - Indicates whether or not the visible area is
    ///                           sampled or not (default to `true`).
    pub fn new(alpha_x: Float, alpha_y: Float, sample_visible_area: bool) -> Self {
        Self {
            sample_visible_area,
            alpha_x: max(0.001, alpha_x),
            alpha_y: max(0.001, alpha_y),
        }
    }

    /// Maps scalar roughness parameter in [0, 1] to alpha values where
    /// values close to 0 are near-perfect specular reflection.
    ///
    /// * `roughness` - Roughness parameter value.
    pub fn roughness_to_alpha(roughness: Float) -> Float {
        let roughness = max(roughness, 1e-3);
        let x = roughness.ln();
        1.62142
            + 0.819955 * x
            + 0.1734 * x * x
            + 0.0171201 * x * x * x
            + 0.000640711 * x * x * x * x
    }
}

impl MicrofacetDistribution for TrowbridgeReitzDistribution {
    /// Returns whether or not the visible area is sampled or not.
    fn get_sample_visible_area(&self) -> bool {
        self.sample_visible_area
    }

    /// Return the differential area of microfacets oriented with the surface
    /// normal `wh`.
    ///
    /// * `wh` - A sample normal from the distrubition of normal vectors.
    #[rustfmt::skip]
    fn d(&self, wh: &Vector3f) -> Float {
        let tan2_theta = tan_2_theta(wh);
        if tan2_theta.is_infinite() {
            0.0
        } else {
            let cos4_theta = cos_2_theta(wh) * cos_2_theta(wh);
            let e =
                (cos_2_phi(wh) / (self.alpha_x * self.alpha_x) + 
                 sin_2_phi(wh) / (self.alpha_y * self.alpha_y)) *
                tan2_theta;
            1.0 / (PI * self.alpha_x * self.alpha_y * cos4_theta * (1.0 + e) * (1.0 + e))
        }
    }

    /// Returns the invisible masked microfacet area per visible microfacet area.
    ///
    /// * `w` - The direction from camera/viewer.
    #[rustfmt::skip]
    fn lambda(&self, w: &Vector3f) -> Float {
        let abs_tan_theta = abs(tan_theta(w));
        if abs_tan_theta.is_infinite() { 
            0.0 
        } else {
            // Compute _alpha_ for direction _w_
            let alpha = (cos_2_phi(w) * self.alpha_x * self.alpha_x + 
                         sin_2_phi(w) * self.alpha_y * self.alpha_y).sqrt();
            let alpha2_tan2_theta = (alpha * abs_tan_theta) * (alpha * abs_tan_theta);
            (-1.0 + (1.0 + alpha2_tan2_theta).sqrt()) / 2.0
        }
    }

    /// Returns a sample from the distribution of normal vectors.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    #[rustfmt::skip]
    fn sample_wh(&self, wo: &Vector3f, u: &Point2f) -> Vector3f {
        if !self.sample_visible_area {
            let mut phi = TWO_PI * u[1];
            let cos_theta = if self.alpha_x == self.alpha_y {
                let tan_theta2 = self.alpha_x * self.alpha_x * u[0] / (1.0 - u[0]);
                1.0 / (1.0 + tan_theta2).sqrt()
            } else {
                phi = atan(self.alpha_y / self.alpha_x * tan(TWO_PI * u[1] + 0.5 * PI));
                if u[1] > 0.5 {
                    phi += PI;
                }
                let sin_phi = sin(phi);
                let cos_phi = cos(phi);
                let alphax2 = self.alpha_x * self.alpha_x;
                let alphay2 = self.alpha_y * self.alpha_y;
                let alpha2 = 1.0 / (cos_phi * cos_phi / alphax2 + sin_phi * sin_phi / alphay2);
                let tan_theta2 = alpha2 * u[0] / (1.0 - u[0]);
                1.0 / (1.0 + tan_theta2).sqrt()
            };
            let sin_theta = max(0.0, 1.0 - cos_theta * cos_theta).sqrt();
            let wh = spherical_direction(sin_theta, cos_theta, phi);
            if !same_hemisphere(wo, &wh) { -wh } else { wh }
        } else {
            let flip = wo.z < 0.0;
            let wo_copy = if flip { -(*wo) } else { *wo };
            let wh = trowbridge_reitz_sample(&wo_copy, self.alpha_x, self.alpha_y, u[0], u[1]);
            if flip { -wh } else { wh }
        }
    }
}

/// Helper function for sampling visible area of normals.
///
/// * `cos_theta` - Cosine of the angle θ measured from the incident direction
///                 to the z-axis.
/// * `u1`        - The uniform random value.
/// * `u2`        - The uniform random value.
fn trowbridge_reitz_sample_11(cos_theta: Float, u1: Float, u2: Float) -> (Float, Float) {
    // special case (normal incidence)
    if cos_theta > 0.9999 {
        let r = (u1 / (1.0 - u1)).sqrt();
        let phi = 6.28318530718 * u2; // TODO: Why not use TWO_PI * u2.
        let slope_x = r * cos(phi);
        let slope_y = r * sin(phi);
        return (slope_x, slope_y);
    }

    let sin_theta = (max(0.0, 1.0 - cos_theta * cos_theta)).sqrt();
    let tan_theta = sin_theta / cos_theta;
    let a = 1.0 / tan_theta;
    let g1 = 2.0 / (1.0 + (1.0 + 1.0 / (a * a)).sqrt());

    // Sample slope_x.
    let a = 2.0 * u1 / g1 - 1.0;
    let mut tmp = 1.0 / (a * a - 1.0);
    if tmp > 1e10 {
        tmp = 1e10;
    }

    let b = tan_theta;
    let d = max(b * b * tmp * tmp - (a * a - b * b) * tmp, 0.0).sqrt();
    let slope_x_1 = b * tmp - d;
    let slope_x_2 = b * tmp + d;
    let slope_x = if a < 0.0 || slope_x_2 > 1.0 / tan_theta { 
        slope_x_1 
    } else {
        slope_x_2
    };

    // sample slope_y
    let (s, u2) = if u2 > 0.5 {
        (1.0, 2.0 * (u2 - 0.5))
    } else {
        (-1.0, 2.0 * (0.5 - u2))
    };
    let z =
        (u2 * (u2 * (u2 * 0.27385 - 0.73369) + 0.46341)) /
        (u2 * (u2 * (u2 * 0.093073 + 0.309420) - 1.000000) + 0.597999);
    let slope_y = s * z * (1.0 + slope_x * slope_x).sqrt();

    assert!(slope_y.is_finite());
    assert!(!slope_y.is_nan());

    (slope_x, slope_y)
}

/// Helper function for sampling visible area of normals.
///
/// * `wi`      - Incident direction.
/// * `alpha_x` - For microfacets oriented perpendicular to the x-axis and where
///               α = sqrt(2) * σ and σ is the RMS slope of microfacets.
/// * `alpha_y` - For microfacets oriented perpendicular to the y-axis and where
///               α = sqrt(2) * σ and σ is the RMS slope of microfacets.
/// * `u1`      - The uniform random value.
/// * `u2`      - The uniform random value.
fn trowbridge_reitz_sample(wi: &Vector3f, alpha_x: Float, alpha_y: Float, u1: Float, u2: Float) -> Vector3f {
    // 1. Stretch wi.
    let wi_stretched = Vector3f::new(alpha_x * wi.x, alpha_y * wi.y, wi.z).normalize();

    // 2. Simulate P22_{wi}(x_slope, y_slope, 1, 1)
    let (mut slope_x, mut slope_y) = trowbridge_reitz_sample_11(cos_theta(&wi_stretched), u1, u2);

    // 3. Rotate.
    let tmp = cos_phi(&wi_stretched) * slope_x - sin_phi(&wi_stretched) * slope_y;
    slope_y = sin_phi(&wi_stretched) * slope_x + cos_phi(&wi_stretched) * slope_y;
    slope_x = tmp;

    // 4. Unstretch.
    slope_x = alpha_x * slope_x;
    slope_y = alpha_y * slope_y;

    // 5. Compute normal.
    Vector3f::new(-slope_x, -slope_y, 1.0).normalize()
}