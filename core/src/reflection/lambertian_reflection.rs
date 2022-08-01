//! Lambertian Reflection

#![allow(dead_code)]

use super::*;
use std::fmt;

/// BRDF for the Lambertian model for perfect diffuse surfaces that scatters
/// incident illumination equally in all directions.
#[derive(Clone)]
pub struct LambertianReflection {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// Reflectance spectrum which gives the fraction of incident light that
    /// is scattered.
    r: Spectrum,
}

impl LambertianReflection {
    /// Creates a new instance of `LambertianReflection`.
    ///
    /// * `r`     - Reflectance spectrum which gives the fraction of incident
    ///             light that is scattered.
    pub fn new(r: Spectrum) -> BxDF {
        let model = Self {
            bxdf_type: BxDFType::BSDF_REFLECTION | BxDFType::BSDF_DIFFUSE,
            r,
        };
        BxDF::LambertianReflection(model)
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
    pub fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        self.r * INV_PI
    }

    /// Computes the hemispherical-directional reflectance function ρ.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - Samples used by Monte Carlo algorithm.
    pub fn rho_hd(&self, _wo: &Vector3f, _u: &[Point2f]) -> Spectrum {
        self.r
    }

    /// Computes the hemispherical-hemispherical-directional reflectance function ρ.
    ///
    /// * `u1` - Samples used b Monte Carlo algorithm.
    /// * `u2` - Samples used b Monte Carlo algorithm.
    pub fn rho_hh(&self, u1: &[Point2f], u2: &[Point2f]) -> Spectrum {
        assert!(u1.len() == u2.len());
        self.r
    }
}

impl fmt::Display for LambertianReflection {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "LambertianReflection {{ bxdf_type: {}, r: {} }}",
            self.bxdf_type, self.r
        )
    }
}
