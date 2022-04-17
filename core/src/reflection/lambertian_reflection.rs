//! Lambertian Reflection

#![allow(dead_code)]

use super::*;
use bumpalo::Bump;
use std::fmt;

/// BRDF for the Lambertian model for perfect diffuse surfaces that scatters
/// incident illumination equally in all directions.
pub struct LambertianReflection {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// Reflectance spectrum which gives the fraction of incident light that
    /// is scattered.
    r: Spectrum,
}

impl LambertianReflection {
    /// Allocate a new instance of `LambertianReflection`.
    ///
    /// * `arena` - The arena for memory allocations.
    /// * `r`     - Reflectance spectrum which gives the fraction of incident
    ///             light that is scattered.
    #[allow(clippy::mut_from_ref)]
    pub fn alloc<'arena>(arena: &'arena Bump, r: Spectrum) -> &'arena mut BxDF {
        let model = arena.alloc(Self {
            bxdf_type: BxDFType::BSDF_REFLECTION | BxDFType::BSDF_DIFFUSE,
            r,
        });
        arena.alloc(BxDF::LambertianReflection(model))
    }

    /// Clone into a newly allocated a new instance of `LambertianReflection`.
    ///
    /// * `arena` - The arena for memory allocations.
    #[allow(clippy::mut_from_ref)]
    pub fn clone_alloc<'arena>(&self, arena: &'arena Bump) -> &'arena mut BxDF<'arena> {
        let model = arena.alloc(Self {
            bxdf_type: self.bxdf_type,
            r: self.r,
        });
        arena.alloc(BxDF::LambertianReflection(model))
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
