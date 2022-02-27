//! KdSubsurface Material

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

/// Implements subsurface scattering material that allows the scattering properties
/// to vary as a function of the position on the surface.
pub struct KdSubsurfaceMaterial {
    /// Scale factor for absorption and scattering coefficients.
    scale: Float,

    /// Coefficient of diffuse reflection.
    kd: ArcTexture<Spectrum>,

    /// Coefficient of glossy reflection.
    kr: ArcTexture<Spectrum>,

    /// Coefficient of glossy transmission.
    kt: ArcTexture<Spectrum>,

    /// Mean free path `1/σt`.
    mfp: ArcTexture<Spectrum>,

    /// Optional microfacet roughness in the u direction. If zero, perfect
    /// specular reflection is modeled.
    u_roughness: ArcTexture<Float>,

    /// Optional microfacet roughness in the v direction. If zero, perfect
    /// specular reflection is modeled.
    v_roughness: ArcTexture<Float>,

    /// If true, roughness values are expected to be in the range [0,1], and are
    /// remapped to microfacet distribution function parameter values that range
    /// from near-perfect-specular at 0 to very rough at 1. Otherwise the
    /// roughness parameters are used directly for the alpha parameters of the
    /// microfacet distribution function.
    remap_roughness: bool,

    /// Index of refraction of the scattering medium.
    eta: Float,

    /// Table for scattering profile data.
    table: Arc<BSSRDFTable>,

    /// Bump map.
    bump_map: Option<ArcTexture<Float>>,
}

impl KdSubsurfaceMaterial {
    /// Create a new `SubsurfaceMaterial`.
    ///
    /// * `scale`           - Scale factor for absorption and scattering coefficients.
    /// * `kd`              - Coefficient of diffuse reflection.
    /// * `kr`              - Coefficient of glossy reflection.
    /// * `kt`              - Coefficient of glossy transmission.
    /// * `mfp`             - Mean free path `1/σt`.
    /// * `g`               - Asymmetry parameter for Henyey-Greenstein phase function.
    /// * `eta`             - Index of refraction of the scattering medium.
    /// * `u_roughness`     - Optional microfacet roughness in the u direction.
    ///                       If zero, perfect specular reflection is modeled.
    /// * `v_roughness`     - Optional microfacet roughness in the v direction.
    ///                       If zero, perfect specular reflection is modeled.
    /// * `remap_roughness` - If true, roughness values are expected to be in
    ///                       the range [0,1], and are remapped to microfacet
    ///                       distribution function parameter values that range
    ///                       from near-perfect-specular at 0 to very rough at 1.
    ///                       Otherwise the roughness parameters are used directly
    ///                       for the alpha parameters of the microfacet
    ///                       distribution function.
    /// * `bump_map`        - Optional bump map.
    pub fn new(
        scale: Float,
        kd: ArcTexture<Spectrum>,
        kr: ArcTexture<Spectrum>,
        kt: ArcTexture<Spectrum>,
        mfp: ArcTexture<Spectrum>,
        g: Float,
        eta: Float,
        u_roughness: ArcTexture<Float>,
        v_roughness: ArcTexture<Float>,
        remap_roughness: bool,
        bump_map: Option<ArcTexture<Float>>,
    ) -> Self {
        let mut table = BSSRDFTable::new(100, 64);
        table.compute_beam_diffusion(g, eta);

        Self {
            scale,
            kd,
            kr,
            kt,
            mfp,
            eta,
            table: Arc::new(table),
            u_roughness,
            v_roughness,
            remap_roughness,
            bump_map,
        }
    }
}

impl Material for KdSubsurfaceMaterial {
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
        allow_multiple_lobes: bool,
        bsdf: &mut Option<&'arena mut BSDF<'scene>>,
        bssrdf: &mut Option<BSSRDF>,
    ) where
        'arena: 'scene,
    {
        // Perform bump mapping with `bump_map`, if present.
        if let Some(bump_map) = &self.bump_map {
            Material::bump(self, bump_map, si);
        }

        // Initialize BSDF for _SubsurfaceMaterial_
        let r = self.kr.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        let t = self.kt.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        let mut urough = self.u_roughness.evaluate(&si.hit, &si.uv, &si.der);
        let mut vrough = self.v_roughness.evaluate(&si.hit, &si.uv, &si.der);

        let result = BSDF::alloc(arena, &si.hit, &si.shading, Some(self.eta));

        // Evaluate textures for `SubsurfaceMaterial` material and allocate BRDF.
        if r.is_black() && t.is_black() {
            return;
        }

        let is_specular = urough == 0.0 && vrough == 0.0;
        if is_specular && allow_multiple_lobes {
            result.add(FresnelSpecular::alloc(arena, r, t, 1.0, self.eta, mode));
        } else {
            if self.remap_roughness {
                urough = TrowbridgeReitzDistribution::roughness_to_alpha(urough);
                vrough = TrowbridgeReitzDistribution::roughness_to_alpha(vrough);
            }

            if is_specular {
                if !r.is_black() {
                    let fresnel = FresnelDielectric::alloc(arena, 1.0, self.eta);
                    result.add(SpecularReflection::alloc(arena, r, fresnel));
                }
                if !t.is_black() {
                    result.add(SpecularTransmission::alloc(arena, t, 1.0, self.eta, mode));
                }
            } else {
                if !r.is_black() {
                    let fresnel = FresnelDielectric::alloc(arena, 1.0, self.eta);
                    let distrib = TrowbridgeReitzDistribution::alloc(arena, urough, vrough, true);
                    result.add(MicrofacetReflection::alloc(arena, r, distrib, fresnel));
                }
                if !t.is_black() {
                    let distrib = TrowbridgeReitzDistribution::alloc(arena, urough, vrough, true);
                    result.add(MicrofacetTransmission::alloc(
                        arena, t, distrib, 1.0, self.eta, mode,
                    ));
                }
            }
        }

        let mfree = self.scale * self.mfp.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        let kd = self.kd.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        let (sig_a, sig_s) = self.table.subsurface_from_diffuse(&kd, &mfree);

        *bsdf = Some(result);
        *bssrdf = Some(BSSRDF::Tabulated {
            eta: self.eta,
            sigma_a: sig_a,
            sigma_s: sig_s,
            table: Arc::clone(&self.table),
        });
    }
}

impl From<&TextureParams> for KdSubsurfaceMaterial {
    /// Create a Subsurface material from given parameter set.
    ///
    /// * `tp` - Texture parameter set.
    fn from(tp: &TextureParams) -> Self {
        let eta = tp.find_float("eta", 1.33);
        let scale = tp.find_float("scale", 1.0);
        let g = tp.find_float("g", 0.0);

        let kd = tp.get_spectrum_texture_or_else(
            "Kd",
            RGBSpectrum::from_rgb(&[0.5, 0.5, 0.5], None),
            |v| Arc::new(ConstantTexture::new(v)),
        );
        let kr = tp.get_spectrum_texture_or_else("Kr", Spectrum::ONE, |v| {
            Arc::new(ConstantTexture::new(v))
        });
        let kt = tp.get_spectrum_texture_or_else("Kt", Spectrum::ONE, |v| {
            Arc::new(ConstantTexture::new(v))
        });
        let mfp = tp.get_spectrum_texture_or_else("mfp", Spectrum::ONE, |v| {
            Arc::new(ConstantTexture::new(v))
        });

        let u_roughness =
            tp.get_float_texture_or_else("uroughness", 0.0, |v| Arc::new(ConstantTexture::new(v)));

        let v_roughness =
            tp.get_float_texture_or_else("vroughness", 0.0, |v| Arc::new(ConstantTexture::new(v)));

        let bump_map = tp.get_float_texture_or_none("bumpmap");
        let remap_roughness = tp.find_bool("remaproughness", true);

        KdSubsurfaceMaterial::new(
            scale,
            kd,
            kr,
            kt,
            mfp,
            g,
            eta,
            u_roughness,
            v_roughness,
            remap_roughness,
            bump_map,
        )
    }
}
