//! Plastic Material

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

/// Implements plastic material as a mixture of a diffuse and glossy scattering
/// function with parameters controlling the particular colors and specular
/// highlight size.
pub struct PlasticMaterial {
    /// Spectral diffuse reflection.
    kd: ArcTexture<Spectrum>,

    /// Spectral specular reflection.
    ks: ArcTexture<Spectrum>,

    /// Roughness.
    roughness: ArcTexture<Float>,

    /// Bump map.
    bump_map: Option<ArcTexture<Float>>,

    /// Remap roughness value to [0, 1] where higher values represent larger
    /// highlights. If this is `false`, use the microfacet distributions `alpha`
    /// parameter.
    remap_roughness: bool,
}

impl PlasticMaterial {
    /// Create a new `PlasticMaterial`.
    ///
    /// * `kd`              - Spectral diffuse reflection.
    /// * `ks`              - Spectral specular reflection.
    /// * `roughness`       - Roughness.
    /// * `remap_roughness` - Remap roughness value to [0, 1] where higher values
    ///                       represent larger highlights. If this is `false`,
    ///                       use the microfacet distributions `alpha` parameter.
    /// * `bump_map`        - Optional bump map.
    pub fn new(
        kd: ArcTexture<Spectrum>,
        ks: ArcTexture<Spectrum>,
        roughness: ArcTexture<Float>,
        remap_roughness: bool,
        bump_map: Option<ArcTexture<Float>>,
    ) -> Self {
        Self {
            kd,
            ks,
            roughness,
            remap_roughness,
            bump_map,
        }
    }
}

impl Material for PlasticMaterial {
    /// Initializes representations of the light-scattering properties of the
    /// material at the intersection point on the surface.
    ///
    /// * `si`                   - The surface interaction at the intersection.
    /// * `mode`                 - Transport mode (ignored).
    /// * `allow_multiple_lobes` - Indicates whether the material should use
    ///                            BxDFs that aggregate multiple types of
    ///                            scattering into a single BxDF when such BxDFs
    ///                            are available (ignored).
    /// * `bsdf`                 - The computed BSDF.
    /// * `bssrdf`               - The computed BSSSRDF.
    fn compute_scattering_functions<'scene>(
        &self,
        si: &mut SurfaceInteraction<'scene>,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
        bsdf: &mut Option<BSDF>,
        bssrdf: &mut Option<BSSRDF>,
    ) {
        // Perform bump mapping with `bump_map`, if present.
        if let Some(bump_map) = &self.bump_map {
            Material::bump(self, bump_map, si);
        }

        let mut result = BSDF::new(&si.hit, &si.shading, None);

        // Initialize diffuse component of plastic material.
        let kd = self.kd.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        if !kd.is_black() {
            result.add(LambertianReflection::new(kd));
        }

        // Initialize specular component of plastic material.
        let ks = self.ks.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        if !ks.is_black() {
            let fresnel = FresnelDielectric::new(1.5, 1.0);

            // Create microfacet distribution for plastic material.
            let mut rough = self.roughness.evaluate(&si.hit, &si.uv, &si.der);
            if self.remap_roughness {
                rough = TrowbridgeReitzDistribution::roughness_to_alpha(rough);
            }
            let distrib = TrowbridgeReitzDistribution::new(rough, rough, true);
            result.add(MicrofacetReflection::new(ks, distrib, fresnel));
        }

        *bsdf = Some(result);
        *bssrdf = None;
    }
}

impl From<&TextureParams> for PlasticMaterial {
    /// Create a plastic material from given parameter set.
    ///
    /// * `tp` - Texture parameter set.
    fn from(tp: &TextureParams) -> Self {
        let kd = tp.get_spectrum_texture_or_else("Kd", Spectrum::new(0.25), |v| {
            Arc::new(ConstantTexture::new(v))
        });

        let ks = tp.get_spectrum_texture_or_else("Ks", Spectrum::new(0.25), |v| {
            Arc::new(ConstantTexture::new(v))
        });

        let roughness =
            tp.get_float_texture_or_else("roughness", 0.1, |v| Arc::new(ConstantTexture::new(v)));

        let bump_map = tp.get_float_texture_or_none("bumpmap");
        let remap_roughness = tp.find_bool("remaproughness", true);

        Self::new(kd, ks, roughness, remap_roughness, bump_map)
    }
}
