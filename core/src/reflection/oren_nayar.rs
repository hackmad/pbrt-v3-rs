//! Oren-Nayar Microfacet Model

use super::*;
use std::fmt;

/// BRDF for the Oren-Nayar model for modeling rough surfaces using a microfacet model.
#[derive(Clone)]
pub struct OrenNayar {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// Reflectance spectrum which gives the fraction of incident light that is scattered.
    r: Spectrum,

    /// Precomputed parameter `A` where:
    ///               σ^2
    /// A = 1 - ---------------
    ///          2(σ^2 + 0.33)
    ///
    /// and σ is the Gaussian distribution parameter, the standard deviation of the microfacet orientation angle.
    a: Float,

    /// Precomputed parameter `B` where:
    ///      0.45 * σ^2
    /// B = ------------
    ///      σ^2 + 0.09
    b: Float,
}

impl OrenNayar {
    /// Creates a new instance of `OrenNayar`.
    ///
    /// * `r`     - Reflectance spectrum which gives the fraction of incident light that is scattered.
    /// * `sigma` - The Gaussian distribution parameter, the standard deviation of the microfacet orientation angle (in
    ///             degrees).
    pub fn new(r: Spectrum, sigma: Float) -> BxDF {
        let sigma = sigma.to_radians();
        let sigma2 = sigma * sigma;
        let model = Self {
            bxdf_type: BxDFType::BSDF_REFLECTION | BxDFType::BSDF_DIFFUSE,
            r,
            a: 1.0 - (sigma2 / (2.0 * (sigma2 + 0.33))),
            b: 0.45 * sigma2 / (sigma2 + 0.09),
        };
        BxDF::OrenNayar(model)
    }

    /// Returns the BxDF type.
    pub fn get_type(&self) -> BxDFType {
        self.bxdf_type
    }

    /// Returns the value of the distribution function for the given pair of directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let sin_theta_i = sin_theta(wi);
        let sin_theta_o = sin_theta(wo);

        // Compute cosine term of Oren-Nayar model.
        let max_cos = if (sin_theta_i > 1e-4) && (sin_theta_o > 1e-4) {
            let sin_phi_i = sin_phi(wi);
            let cos_phi_i = cos_phi(wi);
            let sin_phi_o = sin_phi(wo);
            let cos_phi_o = cos_phi(wo);
            let d_cos = cos_phi_i * cos_phi_o + sin_phi_i * sin_phi_o;
            max(0.0, d_cos)
        } else {
            0.0
        };

        // Compute sine and tangent terms of Oren-Nayar model.
        let abs_cos_theta_wo = abs_cos_theta(wo);
        let abs_cos_theta_wi = abs_cos_theta(wi);
        let (sin_alpha, tan_beta) = if abs_cos_theta_wi > abs_cos_theta_wo {
            (sin_theta_o, sin_theta_i / abs_cos_theta_wi)
        } else {
            (sin_theta_i, sin_theta_o / abs_cos_theta_wo)
        };

        self.r * INV_PI * (self.a + self.b * max_cos * sin_alpha * tan_beta)
    }
}

impl fmt::Display for OrenNayar {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "OrenNayar {{ bxdf_type: {}, r: {}, a: {}, b: {} }}",
            self.bxdf_type, self.r, self.a, self.b,
        )
    }
}
