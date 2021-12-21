//! Specular Reflection

use super::*;
use bumpalo::Bump;

/// BRDF for physically plausible specular reflection using Fresnel interface.
#[derive(Clone)]
pub struct SpecularReflection {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// Fresnel interface for dielectrics and conductors.
    fresnel: ArcFresnel,

    /// Spectrum used to scale the reflected colour.
    r: Spectrum,
}

impl SpecularReflection {
    /// Create a new instance of `SpecularReflection`.
    ///
    /// * `fresnel` - Fresnel interface for dielectrics and conductors.
    /// * `r`       - Spectrum used to scale the reflected colour.
    pub fn new(r: Spectrum, fresnel: ArcFresnel) -> Self {
        Self {
            bxdf_type: BxDFType::BSDF_REFLECTION | BxDFType::BSDF_SPECULAR,
            fresnel: Arc::clone(&fresnel),
            r,
        }
    }

    /// Allocate a new instance of `SpecularReflection`.
    ///
    /// * `allocator` - The allocator.
    /// * `fresnel`   - Fresnel interface for dielectrics and conductors.
    /// * `r`         - Spectrum used to scale the reflected colour.
    pub fn alloc(allocator: &Bump, r: Spectrum, fresnel: ArcFresnel) -> BxDF {
        let model = allocator.alloc(Self::new(r, fresnel)).to_owned();
        allocator.alloc(BxDF::SpecularReflection(model)).to_owned()
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
        Spectrum::new(0.0)
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
