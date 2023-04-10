//! Lambertian Reflection

use super::*;
use std::fmt;

/// BRDF for the Lambertian model for perfect transmissive surfaces that scatters incident illumination equally through
/// a surface in all directions.
#[derive(Clone)]
pub struct LambertianTransmission {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// Transmission spectrum which gives the fraction of incident light that is scattered through the surface.
    t: Spectrum,
}

impl LambertianTransmission {
    /// Creates a new instance of `LambertianTransmission`.
    ///
    /// * `t` - Transmission spectrum which gives the fraction of incident light that is scattered through the surface.
    pub fn new(t: Spectrum) -> BxDF {
        let model = Self {
            bxdf_type: BxDFType::BSDF_TRANSMISSION | BxDFType::BSDF_DIFFUSE,
            t,
        };
        BxDF::LambertianTransmission(model)
    }

    /// Returns the BxDF type.
    pub fn get_type(&self) -> BxDFType {
        self.bxdf_type
    }

    /// Returns the value of the distribution function for the given pair of directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    pub fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        self.t * INV_PI
    }

    /// Returns the value of the BxDF given the outgpoing direction.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    pub fn sample_f(&self, wo: &Vector3f, u: &Point2f) -> BxDFSample {
        // Cosine-sample the hemisphere, flipping the direction if necessary.
        let mut wi = cosine_sample_hemisphere(u);
        if wo.z > 0.0 {
            wi.z *= -1.0;
        }
        let pdf = self.pdf(wo, &wi);
        BxDFSample::new(self.f(wo, &wi), pdf, wi, self.get_type())
    }

    /// Evaluates the PDF for the sampling method. Default is based on the cosine-weighted sampling in `BxDF::sample_f()`
    /// default implementation.
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if !same_hemisphere(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0
        }
    }

    /// Computes the hemispherical-directional reflectance function ρ.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - Samples used by Monte Carlo algorithm.
    pub fn rho_hd(&self, _wo: &Vector3f, _u: &[Point2f]) -> Spectrum {
        self.t
    }

    /// Computes the hemispherical-hemispherical-directional reflectance function ρ.
    ///
    /// * `u1` - Samples used b Monte Carlo algorithm.
    /// * `u2` - Samples used b Monte Carlo algorithm.
    pub fn rho_hh(&self, u1: &[Point2f], u2: &[Point2f]) -> Spectrum {
        assert!(u1.len() == u2.len());
        self.t
    }
}

impl fmt::Display for LambertianTransmission {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "LambertianTransmission {{ bxdf_type: {}, t: {} }}",
            self.bxdf_type, self.t
        )
    }
}
