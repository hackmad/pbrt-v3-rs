//! Torrance-Sparrow Microfacet Reflection Model

use super::*;
use crate::microfacet::*;
use std::fmt;

/// BRDF for modeling metallic surfaces using a microfacet distribution.
#[derive(Clone)]
pub struct MicrofacetReflection {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// Fresnel interface for dielectrics and conductors.
    fresnel: Fresnel,

    /// Reflectance spectrum which gives the fraction of incident light that is scattered.
    r: Spectrum,

    /// The microfacet distribution model.
    distribution: MicrofacetDistribution,
}

impl MicrofacetReflection {
    /// Creates a new instance of `MicrofacetReflection`.
    ///
    /// * `r`            - Reflectance spectrum which gives the fraction of incident light that is scattered.
    /// * `distribution` - Microfacet distribution.
    /// * `fresnel`      - Fresnel interface for dielectrics and conductors.
    pub fn new(r: Spectrum, distribution: MicrofacetDistribution, fresnel: Fresnel) -> BxDF {
        let model = Self {
            bxdf_type: BxDFType::BSDF_REFLECTION | BxDFType::BSDF_GLOSSY,
            r,
            distribution,
            fresnel,
        };
        BxDF::MicrofacetReflection(model)
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
        let cos_theta_o = abs_cos_theta(wo);
        let cos_theta_i = abs_cos_theta(wi);
        let wh = wi + wo;

        // Handle degenerate cases for microfacet reflection.
        if (cos_theta_i == 0.0 || cos_theta_o == 0.0) || (wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0) {
            Spectrum::ZERO
        } else {
            let wh = wh.normalize();

            // For the Fresnel call, make sure that wh is in the same hemisphere as the surface normal, so that TIR is
            // handled correctly.
            let f = self
                .fresnel
                .evaluate(wi.dot(&wh.face_forward(&Vector3f::new(0.0, 0.0, 1.0))));

            self.r * self.distribution.d(&wh) * self.distribution.g(wo, wi) * f / (4.0 * cos_theta_i * cos_theta_o)
        }
    }

    /// Returns the value of the BxDF given the outgpoing direction.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    pub fn sample_f(&self, wo: &Vector3f, u: &Point2f) -> BxDFSample {
        // Sample microfacet orientation `wh` and reflected direction `wi`.
        if wo.z == 0.0 {
            BxDFSample::from(self.bxdf_type)
        } else {
            let wh = self.distribution.sample_wh(wo, u);
            if wo.dot(&wh) < 0.0 {
                // Should be rare.
                BxDFSample::from(self.bxdf_type)
            } else {
                let wi = reflect(wo, &wh);
                if !same_hemisphere(wo, &wi) {
                    BxDFSample::new(Spectrum::ZERO, 0.0, wi, self.bxdf_type)
                } else {
                    // Compute PDF of `wi` for microfacet reflection.
                    let pdf = self.distribution.pdf(wo, &wh) / (4.0 * wo.dot(&wh));
                    BxDFSample::new(self.f(wo, &wi), pdf, wi, self.bxdf_type)
                }
            }
        }
    }

    /// Evaluates the PDF for the sampling method. Default is based on the cosine-weighted sampling in `BxDF::sample_f()`
    /// default implementation.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if same_hemisphere(wo, wi) {
            let wh = (wo + wi).normalize();
            self.distribution.pdf(wo, &wh) / (4.0 * wo.dot(&wh))
        } else {
            0.0
        }
    }
}

impl fmt::Display for MicrofacetReflection {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "MicrofacetReflection {{ bxdf_type: {}, fresnel: {}, r: {}, distribution: {} }}",
            self.bxdf_type, self.fresnel, self.r, self.distribution
        )
    }
}
