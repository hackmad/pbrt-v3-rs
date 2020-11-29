//! Spectrum

use super::pbrt::*;

mod cie;
mod common;
mod rgb;
mod rgb_spectrum;
mod sampled_spectrum;

// Re-export
pub use cie::*;
pub use common::CoefficientSpectrum;
pub use rgb::*;
pub use rgb_spectrum::*;
pub use sampled_spectrum::*;
