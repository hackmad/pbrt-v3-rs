//! Ashikhmin-Shirley Microfacet Reflection Model using Fresnel Effects.

use super::*;
use crate::microfacet::*;
use crate::rng::*;
use bumpalo::Bump;

/// BRDF for modeling layered surfaces such as wood using Ashikhmin-Shirley model.
#[derive(Clone)]
pub struct FresnelBlend {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// Reflectance spectrum for diffuse scattering.
    rd: Spectrum,

    /// Reflectance spectrum for specular scattering.
    rs: Spectrum,

    /// The microfacet distribution model.
    distribution: ArcMicrofacetDistribution,
}

impl FresnelBlend {
    /// Create a new instance of `FresnelBlend`.
    ///
    /// * `rd`           - Reflectance spectrum for diffuse scattering.
    /// * `rs`           - Reflectance spectrum for specular scattering.
    /// * `distribution` - Microfacet distribution.
    pub fn new(rd: Spectrum, rs: Spectrum, distribution: ArcMicrofacetDistribution) -> Self {
        Self {
            bxdf_type: BxDFType::BSDF_REFLECTION | BxDFType::BSDF_GLOSSY,
            rd,
            rs,
            distribution: Arc::clone(&distribution),
        }
    }

    /// Allocate a new instance of `FresnelBlend`.
    ///
    /// * `allocator`    - The allocator.
    /// * `rd`           - Reflectance spectrum for diffuse scattering.
    /// * `rs`           - Reflectance spectrum for specular scattering.
    /// * `distribution` - Microfacet distribution.
    pub fn alloc(
        allocator: &Bump,
        rd: Spectrum,
        rs: Spectrum,
        distribution: ArcMicrofacetDistribution,
    ) -> Self {
        allocator.alloc(Self::new(rd, rs, distribution)).to_owned()
    }

    /// Returns the Schlick approximation to the Fresnel reflection equations:
    ///
    /// Fr(cosθ) = R + (1 - R)(1 - cosθ)^5
    ///
    /// * `cos_theta` - Angle made by incident direction.
    fn schlick_fresnel(&self, cos_theta: Float) -> Spectrum {
        return self.rs + (Spectrum::new(1.0) - self.rs) * pow5(1.0 - cos_theta);
    }

    /// Returns the BxDF type.
    pub fn get_type(&self) -> BxDFType {
        self.bxdf_type
    }

    /// Returns the value of the distribution function for the given pair of
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let diffuse = (28.0 / (23.0 * PI))
            * self.rd
            * (Spectrum::new(1.0) - self.rs)
            * (1.0 - pow5(1.0 - 0.5 * abs_cos_theta(wi)))
            * (1.0 - pow5(1.0 - 0.5 * abs_cos_theta(wo)));
        let wh = *wi + *wo;
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            Spectrum::new(0.0)
        } else {
            let wh = wh.normalize();
            let specular = self.distribution.d(&wh)
                / (4.0 * wi.abs_dot(&wh) * max(abs_cos_theta(wi), abs_cos_theta(wo)))
                * self.schlick_fresnel(wi.dot(&wh));
            diffuse + specular
        }
    }

    /// Returns the value of the BxDF given the outgpoing direction.
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    pub fn sample_f(&self, wo: &Vector3f, u: &Point2f) -> BxDFSample {
        let mut u = *u; // Make local copy.

        let wi = if u[0] < 0.5 {
            u[0] = min(2.0 * u[0], ONE_MINUS_EPSILON);

            // Cosine-sample the hemisphere, flipping the direction if necessary.
            let mut wi = cosine_sample_hemisphere(&u);
            if wo.z < 0.0 {
                wi.z *= -1.0
            }
            wi
        } else {
            u[0] = min(2.0 * (u[0] - 0.5), ONE_MINUS_EPSILON);

            // Sample microfacet orientation `wh` and reflected direction `wi`.
            let wh = self.distribution.sample_wh(wo, &u);
            let wi = reflect(wo, &wh);

            if !same_hemisphere(wo, &wi) {
                return BxDFSample::new(Spectrum::new(0.0), 0.0, wi, self.bxdf_type);
            }
            wi
        };
        let pdf = self.pdf(wo, &wi);
        BxDFSample::new(self.f(wo, &wi), pdf, wi, self.bxdf_type)
    }

    /// Evaluates the PDF for the sampling method. Default is based on the
    /// cosine-weighted sampling in `BxDF::sample_f()` default implementation.
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if !same_hemisphere(wo, wi) {
            0.0
        } else {
            let wh = (*wo + *wi).normalize();
            let pdf_wh = self.distribution.pdf(wo, &wh);
            0.5 * (abs_cos_theta(wi) * INV_PI + pdf_wh / (4.0 * wo.dot(&wh)))
        }
    }
}
