//! Glass Material

use core::geometry::*;
use core::material::*;
use core::microfacet::*;
use core::paramset::*;
use core::pbrt::*;
use core::reflection::*;
use core::spectrum::*;
use core::texture::*;
use either::*;
use std::sync::Arc;
use textures::*;

/// Implements purely specular surfaces.
pub struct GlassMaterial {
    ///
    kr: ArcTexture<Spectrum>,

    ///
    kt: ArcTexture<Spectrum>,

    ///
    u_roughness: ArcTexture<Float>,

    ///
    v_roughness: ArcTexture<Float>,

    ///
    index: ArcTexture<Float>,

    /// Bump map.
    bump_map: Option<ArcTexture<Float>>,

    ///
    remap_roughness: bool,
}

impl GlassMaterial {
    /// Create a new `GlassMaterial`.
    ///
    ///
    /// * `kr`              -
    /// * `kt`              -
    /// * `u_roughness`     -
    /// * `v_roughness`     -
    /// * `index`           -
    /// * `bump_map`        - Optional bump map.
    /// * `remap_roughness` -
    pub fn new(
        kr: ArcTexture<Spectrum>,
        kt: ArcTexture<Spectrum>,
        u_roughness: ArcTexture<Float>,
        v_roughness: ArcTexture<Float>,
        index: ArcTexture<Float>,
        bump_map: Option<ArcTexture<Float>>,
        remap_roughness: bool,
    ) -> Self {
        Self {
            kr: Arc::clone(&kr),
            kt: Arc::clone(&kt),
            u_roughness: Arc::clone(&u_roughness),
            v_roughness: Arc::clone(&v_roughness),
            index: Arc::clone(&index),
            bump_map: bump_map.clone(),
            remap_roughness,
        }
    }
}

impl Material for GlassMaterial {
    /// Initializes representations of the light-scattering properties of the
    /// material at the intersection point on the surface.
    ///
    /// * `si`                   - The surface interaction at the intersection.
    /// * `mode`                 - Transport mode (ignored).
    /// * `allow_multiple_lobes` - Indicates whether the material should use
    ///                            BxDFs that aggregate multiple types of
    ///                            scattering into a single BxDF when such BxDFs
    ///                            are available (ignored).
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    ) {
        // Perform bump mapping with `bump_map`, if present.
        if let Some(bump_map) = &self.bump_map {
            Material::bump(self, Arc::clone(bump_map), si);
        }

        let eta = self.index.evaluate(si);
        let mut urough = self.u_roughness.evaluate(si);
        let mut vrough = self.v_roughness.evaluate(si);
        let r = self.kr.evaluate(si).clamp_default();
        let t = self.kt.evaluate(si).clamp_default();

        // Initialize bsdf for smooth or rough dielectric.
        let mut bsdf = BSDF::new(&si, Some(eta));

        // Evaluate textures for `GlassMaterial` material and allocate BRDF
        if !(r.is_black() && t.is_black()) {
            let is_specular = urough == 0.0 && vrough == 0.0;
            if is_specular && allow_multiple_lobes {
                bsdf.add(Arc::new(FresnelSpecular::new(r, t, 1.0, eta, mode)));
            } else {
                if self.remap_roughness {
                    urough = TrowbridgeReitzDistribution::roughness_to_alpha(urough);
                    vrough = TrowbridgeReitzDistribution::roughness_to_alpha(vrough);
                }

                if is_specular {
                    if !r.is_black() {
                        let fresnel = Arc::new(FresnelDielectric::new(1.0, eta));
                        bsdf.add(Arc::new(SpecularReflection::new(r, fresnel)));
                    }

                    if !t.is_black() {
                        bsdf.add(Arc::new(SpecularTransmission::new(t, 1.0, eta, mode)));
                    }
                } else {
                    let distrib: ArcMicrofacetDistribution =
                        Arc::new(TrowbridgeReitzDistribution::new(urough, vrough, true));

                    if !r.is_black() {
                        let fresnel = Arc::new(FresnelDielectric::new(1.0, eta));
                        bsdf.add(Arc::new(MicrofacetReflection::new(
                            r,
                            Arc::clone(&distrib),
                            fresnel,
                        )));
                    }

                    if !t.is_black() {
                        bsdf.add(Arc::new(MicrofacetTransmission::new(
                            t,
                            Arc::clone(&distrib),
                            1.0,
                            eta,
                            mode,
                        )));
                    }
                };
            }
        }

        si.bsdf = Some(bsdf);
    }
}

impl From<&TextureParams> for GlassMaterial {
    /// Create a matte material from given parameter set.
    ///
    /// * `tp` - Texture parameter set.
    fn from(tp: &TextureParams) -> Self {
        let kr = match tp
            .get_spectrum_texture_or_else("Kr", Arc::new(ConstantTexture::new(Spectrum::new(1.0))))
        {
            Left(tex) => tex,
            Right(val) => Arc::new(ConstantTexture::new(val)),
        };

        let kt = match tp
            .get_spectrum_texture_or_else("Kt", Arc::new(ConstantTexture::new(Spectrum::new(1.0))))
        {
            Left(tex) => tex,
            Right(val) => Arc::new(ConstantTexture::new(val)),
        };

        let index = match tp.get_float_texture_or_else("index", Arc::new(ConstantTexture::new(1.5)))
        {
            Left(tex) => tex,
            Right(val) => Arc::new(ConstantTexture::new(val)),
        };

        let eta = match tp.get_float_texture_or_else("eta", index) {
            Left(tex) => tex,
            Right(val) => Arc::new(ConstantTexture::new(val)),
        };

        let u_roughness =
            match tp.get_float_texture_or_else("uroughness", Arc::new(ConstantTexture::new(0.0))) {
                Left(tex) => tex,
                Right(val) => Arc::new(ConstantTexture::new(val)),
            };

        let v_roughness =
            match tp.get_float_texture_or_else("vroughness", Arc::new(ConstantTexture::new(0.0))) {
                Left(tex) => tex,
                Right(val) => Arc::new(ConstantTexture::new(val)),
            };

        let bump_map = tp.get_float_texture("bumpmap");
        let remap_roughness = tp.find_bool("remaproughness", true);

        Self::new(
            kr,
            kt,
            u_roughness,
            v_roughness,
            eta,
            bump_map,
            remap_roughness,
        )
    }
}
