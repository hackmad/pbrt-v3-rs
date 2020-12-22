//! Spectrum

use super::pbrt::*;

mod cie;
mod common;
mod rgb;
mod rgb_spectrum;
mod sampled_spectrum;

// Re-export
pub use cie::*;
pub use common::*;
pub use rgb::*;
pub use rgb_spectrum::*;
pub use sampled_spectrum::*;

/// Default to using `RGBSpectrum` for rendering.
#[cfg(not(feature = "sampled-spectrum"))]
pub type Spectrum = RGBSpectrum;

/// Use `SampledSpectrum` for rendering.
#[cfg(feature = "sampled-spectrum")]
pub type Spectrum = SampledSpectrum;
