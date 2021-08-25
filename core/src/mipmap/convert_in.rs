//! Per-texel Conversion

use crate::pbrt::*;
use crate::spectrum::*;
/// Interface to convert texels into type `Tmemory` for MIPMap generation.
pub trait ConvertIn<Tmemory> {
    /// Convert the texel to the type `M` and apply the scale and inverse
    /// gamma correction to texel values.
    ///
    /// * `scale` - Scale for the texel values.
    /// * `gamma` - Do gamma correction for the texel values.
    fn convert_in(&self, scale: Float, gamma: bool) -> Tmemory;
}

// NOTE: We are doing `impl ConvertIn<RGBSpectrum> for RGBSpectrum` because
// images are generally stored as RGB values.
//
// For `SampledSpectrum` this may not be possible. If there is such thing, we
// should really implement for that case separately.
//
// One possibility would be one where we convert RGBSpectrum -> SampledSpectrum.

impl ConvertIn<RGBSpectrum> for RGBSpectrum {
    /// Convert the texel to the type `Spectrum` and apply the scale and
    /// inverse gamma correction to texel values.
    ///
    /// * `scale` - Scale for the texel values.
    /// * `gamma` - Do gamma correction for the texel values.
    fn convert_in(&self, scale: Float, gamma: bool) -> RGBSpectrum {
        let samples: Vec<Float> = self
            .samples()
            .iter()
            .map(|sample| {
                scale
                    * if gamma {
                        inv_gamma_correct(*sample)
                    } else {
                        *sample
                    }
            })
            .collect();
        Self::from(samples)
    }
}

impl ConvertIn<Float> for RGBSpectrum {
    /// Convert the texel to the type `Float` and apply the scale and inverse
    /// gamma correction to texel values.
    ///
    /// * `scale` - Scale for the texel values.
    /// * `gamma` - Do gamma correction for the texel values.
    fn convert_in(&self, scale: Float, gamma: bool) -> Float {
        scale
            * if gamma {
                inv_gamma_correct(self.y())
            } else {
                self.y()
            }
    }
}
