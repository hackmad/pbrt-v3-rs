//! Image Texture

#![allow(dead_code)]
mod convert_in;
mod tex_info;

use crate::core::geometry::*;
use crate::core::image_io::*;
use crate::core::mipmap::*;
use crate::core::pbrt::*;
use crate::core::spectrum::*;
use crate::core::texture::*;
use convert_in::*;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign};
use std::result::Result;
use std::sync::Arc;

/// Stores an image texture with MIPMaps using texels of type `Tmemory`.
#[derive(Clone)]
pub struct ImageTexture<Tmemory> {
    /// 2D mapping.
    mapping: ArcTextureMapping2D,

    /// The mipmaps.
    mipmap: Arc<MIPMap<Tmemory>>,
}

impl<Tmemory> ImageTexture<Tmemory>
where
    Tmemory: Copy
        + Default
        + Mul<Float, Output = Tmemory>
        + MulAssign<Float>
        + Div<Float, Output = Tmemory>
        + DivAssign<Float>
        + Add<Tmemory, Output = Tmemory>
        + AddAssign
        + Clamp<Float>,
    Spectrum: ConvertIn<Tmemory>,
{
    /// Create a new `ImageTexture<M, R>`.
    ///
    /// * `mapping`          - The 2D mapping.
    /// * `path`             - The path to the image file.
    /// * `filtering_method` - Type of filtering to use for mipmaps.
    /// * `wrap_mode`        - Image wrapping convention.
    /// * `scale`            - Scale for the texel values.
    /// * `gamma`            - Do gamma correction for the texel values.
    pub fn new(
        mapping: ArcTextureMapping2D,
        path: &str,
        filtering_method: FilteringMethod,
        wrap_mode: ImageWrap,
        scale: Float,
        gamma: bool,
    ) -> Self {
        let mipmap = match generate_mipmap(path, filtering_method, wrap_mode, scale, gamma) {
            Ok(m) => m,
            Err(err) => panic!("{:}", err),
        };
        Self { mapping, mipmap }
    }
}

// Implement `Texture<Tresult>` for `ImageTexture<Tmemory>` where `Tresult` is
// the output for texture evaluation to fit with the conventions of PBRT v3.

/// Implement `ImageTexture` stored in MIPMaps as `RGBSpectrum` and evaluate to
/// `Spectrum`.
impl Texture<Spectrum> for ImageTexture<RGBSpectrum> {
    /// Evaluate the texture at surface interaction.
    ///
    /// * `si` - Surface interaction.
    fn evaluate(&self, si: &SurfaceInteraction) -> Spectrum {
        // Get the (s, t) mapping for the intersection.
        let TextureMap2DResult {
            p: st,
            dstdx,
            dstdy,
        } = self.mapping.map(si);

        let mem = self.mipmap.lookup(&st, &dstdx, &dstdy);

        // Convert out to `Spectrum`.
        let rgb = mem.to_rgb();
        let mut ret = Spectrum::default();
        ret.from_rgb(&rgb, SpectrumType::Reflectance);
        ret
    }
}

/// Implement `ImageTexture` stored in MIPMaps as `Float` and evaluate to
/// `Float`.
impl Texture<Float> for ImageTexture<Float> {
    /// Evaluate the texture at surface interaction.
    ///
    /// * `si` - Surface interaction.
    fn evaluate(&self, si: &SurfaceInteraction) -> Float {
        // Get the (s, t) mapping for the intersection.
        let TextureMap2DResult {
            p: st,
            dstdx,
            dstdy,
        } = self.mapping.map(si);

        // Convert out to `Float`.
        self.mipmap.lookup(&st, &dstdx, &dstdy)
    }
}

/// Load an image texture from file and build the MIPMap.
///
/// * `path`             - The path to the image file.
/// * `filtering_method` - Type of filtering to use for mipmaps.
/// * `wrap_mode`        - Image wrapping convention.
/// * `scale`            - Scale for the texel values.
/// * `gamma`            - Do gamma correction for the texel values.
fn generate_mipmap<Tmemory>(
    path: &str,
    filtering_method: FilteringMethod,
    wrap_mode: ImageWrap,
    scale: Float,
    gamma: bool,
) -> Result<Arc<MIPMap<Tmemory>>, String>
where
    Tmemory: Copy
        + Default
        + Mul<Float, Output = Tmemory>
        + MulAssign<Float>
        + Div<Float, Output = Tmemory>
        + DivAssign<Float>
        + Add<Tmemory, Output = Tmemory>
        + AddAssign
        + Clamp<Float>,
    Spectrum: ConvertIn<Tmemory>,
{
    // Create `MipMap` for `filename`.
    let RGBImage {
        pixels: mut texels,
        resolution,
    } = match read_image(path) {
        Ok(img) => img,
        Err(err) => return Err(format!("Error reading texture {}, {:}.", path, err)),
    };

    // Flip image in y; texture coordinate space has (0,0) at the lower
    // left corner.
    for y in 0..resolution.y / 2 {
        for x in 0..resolution.x {
            let o1 = y * resolution.x + x;
            let o2 = (resolution.y - 1 - y) * resolution.x + x;
            let tmp = texels[o1];
            texels[o2] = texels[o1];
            texels[o1] = tmp;
        }
    }

    // Convert texels to type M and create MIPMap.
    let converted_texels: Vec<Tmemory> = texels
        .iter()
        .map(|texel| (*texel).convert_in(scale, gamma))
        .collect();

    Ok(Arc::new(MIPMap::new(
        &resolution,
        &converted_texels,
        filtering_method,
        wrap_mode,
    )))
}
