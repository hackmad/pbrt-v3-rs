//! Mix Material

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::material::*;
use crate::core::paramset::*;
use crate::core::pbrt::*;
use crate::core::reflection::*;
use crate::core::spectrum::*;
use crate::core::texture::*;
use crate::textures::*;
use std::sync::Arc;

/// Combines two materials with varying weights.
pub struct MixMaterial {
    /// First material.
    m1: ArcMaterial,

    /// Second material.
    m2: ArcMaterial,

    /// Texture used to blend between `m1` and `m2`.
    scale: ArcTexture<Spectrum>,
}

impl MixMaterial {
    /// Create a new `MixMaterial`.
    ///
    ///
    /// * `m1`    - First material.
    /// * `m2`    - Second material.
    /// * `scale` - Texture used to blend between `m1` and `m2`.
    pub fn new(m1: ArcMaterial, m2: ArcMaterial, scale: ArcTexture<Spectrum>) -> Self {
        Self {
            m1: m1.clone(),
            m2: m2.clone(),
            scale: scale.clone(),
        }
    }
}

impl Material for MixMaterial {
    /// Initializes representations of the light-scattering properties of the
    /// material at the intersection point on the surface.
    ///
    /// * `si`                   - The surface interaction at the intersection.
    /// * `mode`                 - Transport mode.
    /// * `allow_multiple_lobes` - Indicates whether the material should use
    ///                            BxDFs that aggregate multiple types of
    ///                            scattering into a single BxDF when such BxDFs
    ///                            are available.
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    ) {
        // Compute weights and original BxDFs for mix material.
        let s1 = self.scale.evaluate(si).clamp_default();
        let s2 = (Spectrum::new(1.0) - s1).clamp_default();
        let mut si2 = si.clone();

        self.m1
            .compute_scattering_functions(si, mode, allow_multiple_lobes);
        self.m2
            .compute_scattering_functions(&mut si2, mode, allow_multiple_lobes);

        // Initialize `si.bsdf` with weighted mixture of BxDFs.
        let mut bsdf = BSDF::new(&si.clone(), None);

        if let Some(si_bsdf) = si.bsdf.as_mut() {
            let n1 = si_bsdf.num_components(BxDFType::from(BSDF_ALL));
            for i in 0..n1 {
                bsdf.add(Arc::new(ScaledBxDF::new(si_bsdf.bxdfs[i].clone(), s1)));
            }
        }

        if let Some(si2_bsdf) = si2.bsdf.as_mut() {
            let n2 = si2_bsdf.num_components(BxDFType::from(BSDF_ALL));
            for i in 0..n2 {
                bsdf.add(Arc::new(ScaledBxDF::new(si2_bsdf.bxdfs[i].clone(), s2)));
            }
        }

        si.bsdf = Some(Arc::new(bsdf));
    }
}

impl From<(&TextureParams, ArcMaterial, ArcMaterial)> for MixMaterial {
    /// Create a mix material from given parameter set and materials.
    ///
    /// * `props` - Mix material creation properties.
    fn from(props: (&TextureParams, ArcMaterial, ArcMaterial)) -> Self {
        let (tp, mat1, mat2) = props;
        let scale = tp.get_spectrum_texture_or_else(
            "amount",
            Arc::new(ConstantTexture::new(Spectrum::new(0.5))),
        );
        Self::new(mat1.clone(), mat2.clone(), scale)
    }
}
