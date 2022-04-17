//! Specular Reflection

use super::*;
use bumpalo::Bump;
use std::fmt;

/// BRDF for physically plausible specular reflection using Fresnel interface.
pub struct SpecularReflection<'arena> {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// Fresnel interface for dielectrics and conductors.
    fresnel: &'arena mut Fresnel<'arena>,

    /// Spectrum used to scale the reflected colour.
    r: Spectrum,
}

impl<'arena> SpecularReflection<'arena> {
    /// Allocate a new instance of `SpecularReflection`.
    ///
    /// * `arena`   - The arena for memory allocations.
    /// * `fresnel` - Fresnel interface for dielectrics and conductors.
    /// * `r`       - Spectrum used to scale the reflected colour.
    #[allow(clippy::mut_from_ref)]
    pub fn alloc(
        arena: &'arena Bump,
        r: Spectrum,
        fresnel: &'arena mut Fresnel<'arena>,
    ) -> &'arena mut BxDF<'arena> {
        let model = arena.alloc(Self {
            bxdf_type: BxDFType::BSDF_REFLECTION | BxDFType::BSDF_SPECULAR,
            fresnel,
            r,
        });
        arena.alloc(BxDF::SpecularReflection(model))
    }

    /// Clone into a newly allocated a new instance of `SpecularReflection`.
    ///
    /// * `arena` - The arena for memory allocations.
    #[allow(clippy::mut_from_ref)]
    pub fn clone_alloc(&self, arena: &'arena Bump) -> &'arena mut BxDF<'arena> {
        let fresnel = self.fresnel.clone_alloc(arena);
        let model = arena.alloc(Self {
            bxdf_type: self.bxdf_type,
            fresnel,
            r: self.r,
        });
        arena.alloc(BxDF::SpecularReflection(model))
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
        // No scattering is returned.
        Spectrum::ZERO
    }

    /// Returns the value of the BxDF given the outgpoing direction.
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    pub fn sample_f(&self, wo: &Vector3f, _u: &Point2f) -> BxDFSample {
        // Compute perfect specular reflection direction.
        let wi = Vector3f::new(-wo.x, -wo.y, wo.z);
        let pdf = 1.0;
        let s = self.fresnel.evaluate(cos_theta(&wi)) * self.r / abs_cos_theta(&wi);
        BxDFSample::new(s, pdf, wi, self.bxdf_type)
    }

    /// Evaluates the PDF for the sampling method. Default is based on the
    /// cosine-weighted sampling in `BxDF::sample_f()` default implementation.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    pub fn pdf(&self, _wo: &Vector3f, _wi: &Vector3f) -> Float {
        0.0
    }
}

impl<'arena> fmt::Display for SpecularReflection<'arena> {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SpecularReflection {{ bxdf_type: {}, fresnel: {}, r: {} }}",
            self.bxdf_type, self.fresnel, self.r,
        )
    }
}
