//! Lambertian Reflection

#![allow(dead_code)]

use super::*;

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
    /// Create a new instance of `LambertianReflection`.
    ///
    /// * `r` - Reflectance spectrum which gives the fraction of incident light
    ///         that is scattered.
    pub fn new(r: Spectrum) -> Self {
        Self {
            bxdf_type: BSDF_REFLECTION | BSDF_DIFFUSE,
            r,
        }
    }
}

impl BxDF for LambertianReflection {
    /// Returns the BxDF type.
    fn get_type(&self) -> BxDFType {
        self.bxdf_type
    }

    /// Returns the value of the distribution function for the given pair of
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        self.r * INV_PI
    }

    /// Computes the hemispherical-directional reflectance function ρ.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - Samples used by Monte Carlo algorithm.
    fn rho_hd(&self, _wo: &Vector3f, _u: &[Point2f]) -> Spectrum {
        self.r
    }

    /// Computes the hemispherical-hemispherical-directional reflectance function ρ.
    ///
    /// * `u1` - Samples used b Monte Carlo algorithm.
    /// * `u2` - Samples used b Monte Carlo algorithm.
    fn rho_hh(&self, u1: &[Point2f], u2: &[Point2f]) -> Spectrum {
        assert!(u1.len() == u2.len());
        self.r
    }
}
