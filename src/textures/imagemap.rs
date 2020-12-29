//! Image Texture

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::mipmap::*;
use crate::core::pbrt::*;
use crate::core::spectrum::*;
use crate::core::texture::*;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign};
use std::sync::Arc;

/// Stores an image texture with MIPMaps using texels of type `Tmemory`.
#[derive(Clone)]
pub struct ImageTexture<Tmemory>
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
    /// 2D mapping.
    mapping: ArcTextureMapping2D,

    /// The mipmaps.
    mipmap: Arc<MIPMap<Tmemory>>,
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
