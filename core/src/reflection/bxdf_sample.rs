//! BxDF Sample

use super::*;
/// Stores evaluation of BxDF samples.
#[derive(Copy, Clone, Default)]
pub struct BxDFSample {
    /// The sample value.
    pub f: Spectrum,

    /// The value of the PDF.
    pub pdf: Float,

    /// The sampled inbound direction.
    pub wi: Vector3f,

    /// The type of BxDF.
    pub bxdf_type: BxDFType,
}

impl BxDFSample {
    /// Create a new `BxDFSample`.
    ///
    /// * `f`            - The sample value.
    /// * `pdf`          - The value of the PDF.
    /// * `wi`           - The sampled inbound direction.
    /// * `sampled_type` - The type of BxDF.
    pub fn new(f: Spectrum, pdf: Float, wi: Vector3f, bxdf_type: BxDFType) -> Self {
        Self {
            f,
            pdf,
            wi,
            bxdf_type,
        }
    }
}
impl From<BxDFType> for BxDFSample {
    /// Create a `BxDFSample` with just the BxDF type and other fields set to
    /// default.
    ///
    /// * `sampled_type` - The type of BxDF.
    fn from(sampled_type: BxDFType) -> Self {
        Self::new(Spectrum::default(), 0.0, Vector3f::default(), sampled_type)
    }
}
