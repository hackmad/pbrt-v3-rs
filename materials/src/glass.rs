//! Glass Material

use core::bssrdf::*;
use core::interaction::*;
use core::material::*;
use core::microfacet::*;
use core::paramset::*;
use core::pbrt::*;
use core::reflection::*;
use core::spectrum::*;
use core::texture::*;
use std::sync::Arc;
use textures::*;

/// Implements perfect or glossy specular reflection and transmission, weighted by Fresnel terms for accurate
/// angular-dependent variation.
pub struct GlassMaterial {
    /// Reflectivity of the surface.
    kr: ArcTexture<Spectrum>,

    /// Transmissibity of the surface.
    kt: ArcTexture<Spectrum>,

    /// Microfacet roughness in the u direction. If zero, perfect specular
    /// reflection is modeled.
    u_roughness: ArcTexture<Float>,

    /// Microfacet roughness in the v direction. If zero, perfect specular
    /// reflection is modeled.
    v_roughness: ArcTexture<Float>,

    /// The index of refraction of the inside of the object. Implicitly assumes that the exterior of objects is a
    /// vacuum, with IOR of 1.
    index: ArcTexture<Float>,

    /// Bump map.
    bump_map: Option<ArcTexture<Float>>,

    /// If true, roughness values are expected to be in the range [0,1], and are remapped to microfacet distribution
    /// function parameter values that range from near-perfect-specular at 0 to very rough at 1. Otherwise the roughness
    /// parameters are used directly for the alpha parameters of the microfacet distribution function.
    remap_roughness: bool,
}

impl GlassMaterial {
    /// Create a new `GlassMaterial`.
    ///
    /// * `kr`              - Reflectivity of the surface.
    /// * `kt`              - Transmissibity of the surface.
    /// * `u_roughness`     - Microfacet roughness in the u direction. If zero, perfect specular reflection is modeled.
    /// * `v_roughness`     - Microfacet roughness in the u direction. If zero, perfect specular reflection is modeled.
    /// * `index`           - Bump map.
    /// * `bump_map`        - Optional bump map.
    /// * `remap_roughness` - If true, roughness values are expected to be in the range [0,1], and are remapped to
    ///                       microfacet distribution function parameter values that range from near-perfect-specular at
    ///                       0 to very rough at 1. Otherwise the roughness parameters are used directly for the alpha
    ///                       parameters of the microfacet distribution function.
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
            kr,
            kt,
            u_roughness,
            v_roughness,
            index,
            bump_map,
            remap_roughness,
        }
    }
}

impl Material for GlassMaterial {
    /// Initializes representations of the light-scattering properties of the material at the intersection point on the
    /// surface.
    ///
    /// * `si`                   - The surface interaction at the intersection.
    /// * `mode`                 - Transport mode (ignored).
    /// * `allow_multiple_lobes` - Indicates whether the material should use BxDFs that aggregate multiple types of
    ///                            scattering into a single BxDF when such BxDFs are available (ignored).
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
        // Perform bump mapping with `bump_map`, if present.
        if let Some(bump_map) = &self.bump_map {
            Material::bump(self, bump_map, si);
        }

        let eta = self.index.evaluate(&si.hit, &si.uv, &si.der);
        let mut urough = self.u_roughness.evaluate(&si.hit, &si.uv, &si.der);
        let mut vrough = self.v_roughness.evaluate(&si.hit, &si.uv, &si.der);
        let r = self.kr.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        let t = self.kt.evaluate(&si.hit, &si.uv, &si.der).clamp_default();

        // Initialize bsdf for smooth or rough dielectric.
        let mut result = BSDF::new(&si.hit, &si.shading, None);

        // Evaluate textures for `GlassMaterial` material and allocate BRDF.
        if !(r.is_black() && t.is_black()) {
            let is_specular = urough == 0.0 && vrough == 0.0;
            if is_specular && allow_multiple_lobes {
                result.add(FresnelSpecular::new(r, t, 1.0, eta, mode));
            } else {
                if self.remap_roughness {
                    urough = TrowbridgeReitzDistribution::roughness_to_alpha(urough);
                    vrough = TrowbridgeReitzDistribution::roughness_to_alpha(vrough);
                }

                if is_specular {
                    if !r.is_black() {
                        let fresnel = FresnelDielectric::new(1.0, eta);
                        result.add(SpecularReflection::new(r, fresnel));
                    }
                    if !t.is_black() {
                        result.add(SpecularTransmission::new(t, 1.0, eta, mode));
                    }
                } else {
                    if !r.is_black() {
                        let distrib = TrowbridgeReitzDistribution::new(urough, vrough, true);
                        let fresnel = FresnelDielectric::new(1.0, eta);
                        result.add(MicrofacetReflection::new(r, distrib, fresnel));
                    }
                    if !t.is_black() {
                        let distrib = TrowbridgeReitzDistribution::new(urough, vrough, true);
                        result.add(MicrofacetTransmission::new(t, distrib, 1.0, eta, mode));
                    }
                };
            }
        }

        *bsdf = Some(result);
        *bssrdf = None;
    }
}

impl From<&TextureParams> for GlassMaterial {
    /// Create a glass material from given parameter set.
    ///
    /// * `tp` - Texture parameter set.
    fn from(tp: &TextureParams) -> Self {
        let kr = tp.get_spectrum_texture_or_else("Kr", Spectrum::ONE, |v| Arc::new(ConstantTexture::new(v)));

        let kt = tp.get_spectrum_texture_or_else("Kt", Spectrum::ONE, |v| Arc::new(ConstantTexture::new(v)));

        let eta = match tp.get_float_texture("eta") {
            None => tp.get_float_texture_or_else("index", 1.5, |v| Arc::new(ConstantTexture::new(v))),
            Some(tex) => tex,
        };

        let u_roughness = tp.get_float_texture_or_else("uroughness", 0.0, |v| Arc::new(ConstantTexture::new(v)));

        let v_roughness = tp.get_float_texture_or_else("vroughness", 0.0, |v| Arc::new(ConstantTexture::new(v)));

        let bump_map = tp.get_float_texture_or_none("bumpmap");
        let remap_roughness = tp.find_bool("remaproughness", true);

        Self::new(kr, kt, u_roughness, v_roughness, eta, bump_map, remap_roughness)
    }
}
