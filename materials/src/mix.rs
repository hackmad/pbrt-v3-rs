//! Mix Material

use core::bssrdf::*;
use core::interaction::*;
use core::material::*;
use core::paramset::*;
use core::pbrt::*;
use core::reflection::*;
use core::spectrum::*;
use core::texture::*;
use std::sync::Arc;
use textures::*;

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
    /// * `m1`    - First material.
    /// * `m2`    - Second material.
    /// * `scale` - Texture used to blend between `m1` and `m2`.
    pub fn new(m1: ArcMaterial, m2: ArcMaterial, scale: ArcTexture<Spectrum>) -> Self {
        Self {
            m1: Arc::clone(&m1),
            m2: Arc::clone(&m2),
            scale: Arc::clone(&scale),
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
    /// * `bsdf`                 - The computed BSDF.
    /// * `bssrdf`               - The computed BSSSRDF.
    fn compute_scattering_functions<'scene>(
        &self,
        si: &mut SurfaceInteraction<'scene>,
        mode: TransportMode,
        allow_multiple_lobes: bool,
        bsdf: &mut Option<BSDF>,
        bssrdf: &mut Option<BSSRDF>,
    ) {
        // Compute weights and original BxDFs for mix material.
        let s1 = self
            .scale
            .evaluate(&si.hit, &si.uv, &si.der)
            .clamp_default();
        let s2 = (Spectrum::ONE - s1).clamp_default();

        let mut si2 = si.clone();

        let mut si_bsdf: Option<BSDF> = None;
        let mut si_bssrdf: Option<BSSRDF> = None;
        self.m1.compute_scattering_functions(
            si,
            mode,
            allow_multiple_lobes,
            &mut si_bsdf,
            &mut si_bssrdf,
        );

        let mut si2_bsdf: Option<BSDF> = None;
        let mut si2_bssrdf: Option<BSSRDF> = None;
        self.m2.compute_scattering_functions(
            &mut si2,
            mode,
            allow_multiple_lobes,
            &mut si2_bsdf,
            &mut si2_bssrdf,
        );

        // Initialize `si.bsdf` with weighted mixture of BxDFs.
        let mut result = BSDF::new(&si.hit, &si.shading, None);

        if let Some(si_bsdf) = si_bsdf {
            let n1 = si_bsdf.num_components(BxDFType::all());
            for i in 0..n1 {
                result.add(ScaledBxDF::new(si_bsdf.bxdfs[i].clone(), s1));
            }
        }

        if let Some(si2_bsdf) = si2_bsdf {
            let n2 = si2_bsdf.num_components(BxDFType::all());
            for i in 0..n2 {
                result.add(ScaledBxDF::new(si2_bsdf.bxdfs[i].clone(), s2));
            }
        }

        *bsdf = Some(result);
        *bssrdf = None; // TODO: Should the BSSRDFs be ignored?
    }
}

impl From<(&TextureParams, ArcMaterial, ArcMaterial)> for MixMaterial {
    /// Create a mix material from given parameter set and materials.
    ///
    /// * `props` - Mix material creation properties.
    fn from(props: (&TextureParams, ArcMaterial, ArcMaterial)) -> Self {
        let (tp, mat1, mat2) = props;

        let scale = tp.get_spectrum_texture_or_else("amount", Spectrum::new(0.5), |v| {
            Arc::new(ConstantTexture::new(v))
        });

        Self::new(Arc::clone(&mat1), Arc::clone(&mat2), scale)
    }
}
