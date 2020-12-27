//! Beckmann–Spizzichino Distribution

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::reflection::*;
use super::MicrofacetDistribution;

/// Implements the Beckmann–Spizzichino distribution which based on Gaussian 
/// distribution of microfacet slopes.
#[derive(Copy, Clone, Default)]
pub struct BeckmannDistribution {
    /// Indicates whether or not the visible area is sampled or not.
    sample_visible_area: bool,

    /// For microfacets oriented perpendicular to the x-axis and where
    /// α = sqrt(2) * σ and σ is the RMS slope of microfacets.
    alpha_x: Float,

    /// For microfacets oriented perpendicular to the y-axis and where
    /// α = sqrt(2) * σ and σ is the RMS slope of microfacets.
    alpha_y: Float,
}

impl BeckmannDistribution {
    /// Create a new `BeckmannDistribution`.
    ///
    /// * `alpha_x`             - For microfacets oriented perpendicular to the
    ///                           x-axis and where α = sqrt(2) * σ and σ is the
    ///                           RMS slope of microfacets.
    /// * `alpha_y`             - For microfacets oriented perpendicular to the
    ///                           y-axis and where α = sqrt(2) * σ and σ is the
    ///                           RMS slope of microfacets.
    /// * `sample_visible_area` - Indicates whether or not the visible area is
    ///                           sampled or not.
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

impl MicrofacetDistribution for BeckmannDistribution {
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
            (-tan2_theta * (cos_2_phi(wh) / (self.alpha_x * self.alpha_x) +
                            sin_2_phi(wh) / (self.alpha_y * self.alpha_y))
            ).exp() /
            (PI * self.alpha_x * self.alpha_y * cos4_theta)
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
            // Compute _alpha_ for direction `w`.
            let alpha = (cos_2_phi(w) * self.alpha_x * self.alpha_x + 
                         sin_2_phi(w) * self.alpha_y * self.alpha_y).sqrt();
            let a = 1.0 / (alpha * abs_tan_theta);
            if a >= 1.6 {
                0.0
            } else {
                (1.0 - 1.259 * a + 0.396 * a * a) / (3.535 * a + 2.181 * a * a)
            }
        }
    }

    /// Returns a sample from the distribution of normal vectors.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    #[rustfmt::skip]
    fn sample_wh(&self, wo: &Vector3f, u: &Point2f) -> Vector3f {
        if !self.sample_visible_area {
             // Sample full distribution of normals for Beckmann distribution.

             // Compute tan^2(theta) and phi for Beckmann distribution sample.
             let (tan2_theta, phi) = if self.alpha_x == self.alpha_y {
                 let log_sample = (1.0 - u[0]).ln();
                 debug_assert!(log_sample.is_finite());

                 let tan2_theta = -self.alpha_x * self.alpha_x * log_sample;
                 let phi = u[1] * TWO_PI;

                 (tan2_theta, phi)
             } else {
                 // Compute tan^2(theta) and phi for anisotropic Beckmann
                 // distribution.
                 let log_sample = (1.0 - u[0]).ln();
                 debug_assert!(log_sample.is_finite());

                 let mut phi = atan(self.alpha_y / self.alpha_x * tan(TWO_PI * u[1] + 0.5 * PI));
                 if u[1] > 0.5 {
                     phi += PI;
                 }
                 let sin_phi = sin(phi);
                 let cos_phi = cos(phi);
                 let alphax2 = self.alpha_x * self.alpha_x;
                 let alphay2 = self.alpha_y * self.alpha_y;
                 let tan2_theta = -log_sample / (cos_phi * cos_phi / alphax2 + sin_phi * sin_phi / alphay2);

                 (tan2_theta, phi)
             };

             // Map sampled Beckmann angles to normal direction _wh_
             let cos_theta = 1.0 / (1.0 + tan2_theta).sqrt();
             let sin_theta = (max(0.0, 1.0 - cos_theta * cos_theta)).sqrt();
             let wh = spherical_direction(sin_theta, cos_theta, phi);
             if !same_hemisphere(wo, &wh) { -wh } else { wh }
         } else {
             // Sample visible area of normals for Beckmann distribution
             let flip = wo.z < 0.0;
             let wo_copy = if flip { -(*wo) } else { *wo };
             let wh = beckmann_sample(&wo_copy, self.alpha_x, self.alpha_y, u[0], u[1]);
             if flip { -wh } else { wh }
         }
    }
}

/// Helper function for sampling visible area of normals.
///
/// * `cos_theta_i` - Cosine of the angle θ measured from the incident direction
///                   to the z-axis.
/// * `u1`          - The uniform random value.
/// * `u2`          - The uniform random value.
fn beckmann_sample_11(cos_theta_i: Float, u1: Float, u2: Float) -> (Float, Float) {
    // Special case (normal incidence).
    if cos_theta_i > 0.9999 {
        let r = (-(1.0 - u1).ln()).sqrt();
        let sin_phi = sin(TWO_PI * u2);
        let cos_phi = cos(TWO_PI * u2);
        let slope_x = r * cos_phi;
        let slope_y = r * sin_phi;
        return (slope_x, slope_y);
    }

    // The original inversion routine from the paper contained discontinuities,
    // which causes issues for QMC integration and techniques like Kelemen-style
    // MLT. The following code performs a numerical inversion with better behavior.
    let sin_theta_i = max(0.0, 1.0 - cos_theta_i * cos_theta_i).sqrt();
    let tan_theta_i = sin_theta_i / cos_theta_i;
    let cot_theta_i = 1.0 / tan_theta_i;

    // Search interval -- everything is parameterized in the Erf() domain.
    let mut a = -1.0;
    let mut c = erf(cot_theta_i);
    let sample_x = max(u1, 1e-6);

    // Start with a good initial guess.
    // Float b = (1-sample_x) * a + sample_x * c;

    // We can do better (inverse of an approximation computed in Mathematica).
    let theta_i = acos(cos_theta_i);
    let fit = 1.0 + theta_i * (-0.876 + theta_i * (0.4265 - 0.0594 * theta_i));
    let mut b = c - (1.0 + c) * (1.0 - sample_x).powf(fit);

    // Normalization factor for the CDF.
    let inv_sqrt_pi = 1.0 / PI.sqrt();
    let normalization =
        1.0 / (1.0 + c + inv_sqrt_pi * tan_theta_i * (-cot_theta_i * cot_theta_i).exp());

    let mut it = 0;
    loop {
        it += 1;

        if it >= 10 {
            break;
        }

        // Bisection criterion -- the oddly-looking Boolean expression are
        // intentional to check for NaNs at little additional cost.
        if !(b >= a && b <= c) {
            b = 0.5 * (a + c);
        }

        // Evaluate the CDF and its derivative (i.e. the density function).
        let inv_erf = erf_inv(b);
        let value = normalization
            * (1.0 + b + inv_sqrt_pi * tan_theta_i * (-inv_erf * inv_erf).exp())
            - sample_x;
        let derivative = normalization * (1.0 - inv_erf * tan_theta_i);

        if abs(value) < 1e-5 {
            break;
        }

        // Update bisection intervals.
        if value > 0.0 {
            c = b;
        } else {
            a = b;
        }

        b -= value / derivative;
    }

    // Now convert back into a slope value.
    let slope_x = erf_inv(b);

    // Simulate Y component.
    let slope_y = erf_inv(2.0 * max(u2, 1e-6) - 1.0);

    assert!(slope_x.is_finite());
    assert!(!slope_x.is_nan());
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
fn beckmann_sample(
    wi: &Vector3f,
    alpha_x: Float,
    alpha_y: Float,
    u1: Float,
    u2: Float,
) -> Vector3f {
    // 1. Stretch wi
    let wi_stretched = Vector3f::new(alpha_x * wi.x, alpha_y * wi.y, wi.z).normalize();

    // 2. Simulate P22_{wi}(x_slope, y_slope, 1, 1).
    let (mut slope_x, mut slope_y) = beckmann_sample_11(cos_theta(&wi_stretched), u1, u2);

    // 3. Rotate
    let tmp = cos_phi(&wi_stretched) * slope_x - sin_phi(&wi_stretched) * slope_y;
    slope_y = sin_phi(&wi_stretched) * slope_x + cos_phi(&wi_stretched) * slope_y;
    slope_x = tmp;

    // 4. unstretch
    slope_x = alpha_x * slope_x;
    slope_y = alpha_y * slope_y;

    // 5. compute normal
    Vector3f::new(-slope_x, -slope_y, 1.0).normalize()
}
