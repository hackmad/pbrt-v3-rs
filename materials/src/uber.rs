//! Uber Material

use bumpalo::Bump;
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

/// Implements a “kitchen sink” material representing the union of many of the
/// preceding materials. This is a highly parameterized material that is
/// particularly useful when converting scenes from other file formats into
/// pbrt’s.
pub struct UberMaterial {
    /// Coefficient of glossy reflection.
    ks: ArcTexture<Spectrum>,

    /// Coefficient of diffuse reflection.
    kd: ArcTexture<Spectrum>,

    /// Coefficient of specular reflection.
    kr: ArcTexture<Spectrum>,

    /// Coefficient of specular transmission.
    kt: ArcTexture<Spectrum>,

    /// The opacity of the surface. Note that when less than one, the uber
    /// material transmits light without refracting it.
    opacity: ArcTexture<Spectrum>,

    /// Roughness of the surface.
    roughness: ArcTexture<Float>,

    /// Optional microfacet roughness in the u direction. If zero, perfect
    /// specular reflection is modeled.
    u_roughness: Option<ArcTexture<Float>>,

    /// Optional microfacet roughness in the v direction. If zero, perfect
    /// specular reflection is modeled.
    v_roughness: Option<ArcTexture<Float>>,

    /// The index of refraction of the inside of the object.
    /// Implicitly assumes that the exterior of objects is a vacuum, with IOR of 1.
    index: ArcTexture<Float>,

    /// If true, roughness values are expected to be in the range [0,1], and are
    /// remapped to microfacet distribution function parameter values that range
    /// from near-perfect-specular at 0 to very rough at 1. Otherwise the
    /// roughness parameters are used directly for the alpha parameters of the
    /// microfacet distribution function.
    remap_roughness: bool,

    /// Bump map.
    bump_map: Option<ArcTexture<Float>>,
}

impl UberMaterial {
    /// Create a new `UberMaterial`.
    ///
    /// * `ks`              - Spectral glossy reflection.
    /// * `kd`              - Spectral diffuse reflection.
    /// * `kr`              - Reflectivity of the surface.
    /// * `kt`              - Transmissibity of the surface.
    /// * `opacity`         - The opacity of the surface. Note that when less
    ///                       than one, the uber material transmits light without
    ///                       refracting it.
    /// * `roughness`       - Optional roughness of the surface.
    /// * `u_roughness`     - Optional microfacet roughness in the u direction.
    ///                       If zero, perfect specular reflection is modeled.
    /// * `v_roughness`     - Optional microfacet roughness in the v direction.
    ///                       If zero, perfect specular reflection is modeled.
    /// * `index`           - The index of refraction of the inside of the object.
    ///                       Implicitly assumes that the exterior of objects is a vacuum, with IOR of 1.
    /// * `remap_roughness` - If true, roughness values are expected to be in
    ///                       the range [0,1], and are remapped to microfacet
    ///                       distribution function parameter values that range
    ///                       from near-perfect-specular at 0 to very rough at 1.
    ///                       Otherwise the roughness parameters are used directly
    ///                       for the alpha parameters of the microfacet
    ///                       distribution function.
    /// * `bump_map`        - Optional bump map.
    pub fn new(
        ks: ArcTexture<Spectrum>,
        kd: ArcTexture<Spectrum>,
        kr: ArcTexture<Spectrum>,
        kt: ArcTexture<Spectrum>,
        opacity: ArcTexture<Spectrum>,
        roughness: ArcTexture<Float>,
        u_roughness: Option<ArcTexture<Float>>,
        v_roughness: Option<ArcTexture<Float>>,
        index: ArcTexture<Float>,
        remap_roughness: bool,
        bump_map: Option<ArcTexture<Float>>,
    ) -> Self {
        Self {
            ks,
            kd,
            kr,
            kt,
            opacity,
            roughness,
            u_roughness,
            v_roughness,
            index,
            remap_roughness,
            bump_map,
        }
    }
}

impl Material for UberMaterial {
    /// Initializes representations of the light-scattering properties of the
    /// material at the intersection point on the surface.
    ///
    /// * `arena`                - The arena for memory allocations.
    /// * `si`                   - The surface interaction at the intersection.
    /// * `mode`                 - Transport mode (ignored).
    /// * `allow_multiple_lobes` - Indicates whether the material should use
    ///                            BxDFs that aggregate multiple types of
    ///                            scattering into a single BxDF when such BxDFs
    ///                            are available (ignored).
    /// * `bsdf`                 - The computed BSDF.
    /// * `bssrdf`               - The computed BSSSRDF.
    fn compute_scattering_functions<'scene, 'arena>(
        &self,
        arena: &'arena Bump,
        si: &mut SurfaceInteraction<'scene>,
        mode: TransportMode,
        _allow_multiple_lobes: bool,
        bsdf: &mut Option<&'arena mut BSDF<'scene>>,
        bssrdf: &mut Option<BSSRDFType>,
    ) where
        'arena: 'scene,
    {
        // Perform bump mapping with `bump_map`, if present.
        if let Some(bump_map) = &self.bump_map {
            Material::bump(self, bump_map, si);
        }

        // Evaluate textures for `UberMaterial` material and allocate BRDF.
        let e = self.index.evaluate(&si.hit, &si.uv, &si.der);

        let op = self
            .opacity
            .evaluate(&si.hit, &si.uv, &si.der)
            .clamp_default();
        let t = (-op + Spectrum::ONE).clamp_default();
        let result = if !t.is_black() {
            let bsdf = BSDF::alloc(arena, &si.hit, &si.shading, Some(1.0));
            let tr = SpecularTransmission::alloc(arena, t, 1.0, 1.0, mode);
            bsdf.add(tr);
            bsdf
        } else {
            BSDF::alloc(arena, &si.hit, &si.shading, Some(e))
        };

        let kd = op * self.kd.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        if !kd.is_black() {
            let diff = LambertianReflection::alloc(arena, kd);
            result.add(diff);
        }

        let ks = op * self.ks.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        if !ks.is_black() {
            let fresnel = FresnelDielectric::alloc(arena, 1.0, e);

            let mut urough = self.u_roughness.as_ref().map_or_else(
                || self.roughness.evaluate(&si.hit, &si.uv, &si.der),
                |roughness| roughness.evaluate(&si.hit, &si.uv, &si.der),
            );

            let mut vrough = self.v_roughness.as_ref().map_or_else(
                || self.roughness.evaluate(&si.hit, &si.uv, &si.der),
                |roughness| roughness.evaluate(&si.hit, &si.uv, &si.der),
            );

            if self.remap_roughness {
                urough = TrowbridgeReitzDistribution::roughness_to_alpha(urough);
                vrough = TrowbridgeReitzDistribution::roughness_to_alpha(vrough);
            }

            let distrib = TrowbridgeReitzDistribution::alloc(arena, urough, vrough, true);
            let spec = MicrofacetReflection::alloc(arena, ks, distrib, fresnel);
            result.add(spec);
        }

        let kr = op * self.kr.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        if !kr.is_black() {
            let fresnel = FresnelDielectric::alloc(arena, 1.0, e);
            result.add(SpecularReflection::alloc(arena, kr, fresnel));
        }

        let kt = op * self.kt.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        if !kt.is_black() {
            result.add(SpecularTransmission::alloc(arena, kt, 1.0, e, mode));
        }

        *bsdf = Some(result);
        *bssrdf = None;
    }
}

impl From<&TextureParams> for UberMaterial {
    /// Create an uber material from given parameter set.
    ///
    /// * `tp` - Texture parameter set.
    fn from(tp: &TextureParams) -> Self {
        let kd = tp.get_spectrum_texture_or_else("Kd", Spectrum::new(0.25), |v| {
            Arc::new(ConstantTexture::new(v))
        });
        let ks = tp.get_spectrum_texture_or_else("Ks", Spectrum::new(0.25), |v| {
            Arc::new(ConstantTexture::new(v))
        });
        let kr = tp.get_spectrum_texture_or_else("Kr", Spectrum::ZERO, |v| {
            Arc::new(ConstantTexture::new(v))
        });
        let kt = tp.get_spectrum_texture_or_else("Kt", Spectrum::ZERO, |v| {
            Arc::new(ConstantTexture::new(v))
        });

        let roughness =
            tp.get_float_texture_or_else("roughness", 0.1, |v| Arc::new(ConstantTexture::new(v)));
        let u_roughness = tp.get_float_texture_or_none("uroughness");
        let v_roughness = tp.get_float_texture_or_none("vroughness");

        let eta = match tp.get_float_texture("eta") {
            None => {
                tp.get_float_texture_or_else("index", 1.5, |v| Arc::new(ConstantTexture::new(v)))
            }
            Some(tex) => tex,
        };

        let opacity = tp.get_spectrum_texture_or_else("opacity", Spectrum::ONE, |v| {
            Arc::new(ConstantTexture::new(v))
        });

        let bump_map = tp.get_float_texture_or_none("bumpmap");
        let remap_roughness = tp.find_bool("remaproughness", true);

        Self::new(
            ks,
            kd,
            kr,
            kt,
            opacity,
            roughness,
            u_roughness,
            v_roughness,
            eta,
            remap_roughness,
            bump_map,
        )
    }
}
