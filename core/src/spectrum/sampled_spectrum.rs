//! Sampled Spectrum.

use super::*;
use crate::pbrt::*;
use std::convert::TryInto;
use std::fmt;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
};

/// Starting wavelength in nm for SPDs.
pub const SAMPLED_LAMBDA_START: usize = 400;

/// Ending wavelength in nm for SPDs.
pub const SAMPLED_LAMBDA_END: usize = 700;

/// Number of spectral samples to use for `SampledSpectrum`.
pub const SPECTRAL_SAMPLES: usize = 60;

/// SampledSpectrum represents an spectral power distribution (SPD) with
/// uniformly spaced samples between a starting and ending wavelength.
///
/// The wavelength range covers human visual range of 400 nm to 700 nm and
/// 60 samples is more than enough to accurately represent complex SPDs.
#[derive(Copy, Clone, PartialEq)]
pub struct SampledSpectrum {
    /// The sampled spectral values.
    c: [Float; SPECTRAL_SAMPLES],
}

impl SampledSpectrum {
    /// Create a new `SampledSpectrum` with a constant value across all
    /// wavelengths.
    ///
    /// * `v` - Constant value.
    pub fn new(v: Float) -> Self {
        let ret = Self {
            c: [v; SPECTRAL_SAMPLES],
        };
        assert!(!ret.has_nans());
        ret
    }

    /// Spectrum with all values set to 0.
    pub const ZERO: Self = Self {
        c: [0.0; SPECTRAL_SAMPLES],
    };

    /// Spectrum with all values set to 1.
    pub const ONE: Self = Self {
        c: [1.0; SPECTRAL_SAMPLES],
    };
}

impl Default for SampledSpectrum {
    /// Return a black `SampledSpectrum`.
    fn default() -> Self {
        Self::ZERO
    }
}

impl From<Float> for SampledSpectrum {
    /// Create a new `SampledSpectrum` with a constant value across all
    /// wavelengths.
    ///
    /// * `v` - Constant value.
    fn from(v: Float) -> Self {
        Self::new(v)
    }
}

impl From<Vec<Float>> for SampledSpectrum {
    /// Create a new `SampledSpectrum` from sampled spectral values.
    ///
    /// * `c` - Sample values.
    fn from(c: Vec<Float>) -> Self {
        let ret = Self {
            c: c.try_into().unwrap_or_else(|v: Vec<Float>| {
                panic!(
                    "Expected a Vec of length {} but it was {}",
                    SPECTRAL_SAMPLES,
                    v.len()
                )
            }),
        };
        assert!(!ret.has_nans());
        ret
    }
}

impl From<&Vec<Sample>> for SampledSpectrum {
    /// Create a new `SampledSpectrum` from sampled spectral samples.
    ///
    /// * `samples` - Samples.
    fn from(samples: &Vec<Sample>) -> Self {
        // Sort samples if unordered.
        let mut sorted_samples = samples.clone();
        if !are_spectrum_samples_sorted(samples) {
            sort_spectrum_samples(&mut sorted_samples);
        };

        let mut c = [0.0; SPECTRAL_SAMPLES];
        for (i, v) in c.iter_mut().enumerate() {
            // Compute average value of given SPD over i^th sample's range.
            let lambda0 = lerp(
                i as Float / SPECTRAL_SAMPLES as Float,
                SAMPLED_LAMBDA_START as Float,
                SAMPLED_LAMBDA_END as Float,
            );
            let lambda1 = lerp(
                (i + 1) as Float / SPECTRAL_SAMPLES as Float,
                SAMPLED_LAMBDA_START as Float,
                SAMPLED_LAMBDA_END as Float,
            );
            *v = average_spectrum_samples(samples, lambda0, lambda1);
        }

        Self { c }
    }
}

impl CoefficientSpectrum for SampledSpectrum {
    /// Returns the stored samples.
    fn samples(&self) -> &[Float] {
        &self.c
    }

    /// Returns stored samples as mutable.
    fn samples_mut(&mut self) -> &mut [Float] {
        &mut self.c
    }

    /// Converts XYZ values to a full SPD.
    ///
    /// * `xyz`           - XYZ colour value.
    /// * `spectrum_type` - Indicates type of colour value. If `None`,
    ///                     defaults to `SpectrumType::Reflectance`.
    fn from_xyz(xyz: &[Float; 3], spectrum_type: Option<SpectrumType>) -> Self {
        Self::from_rgb(&xyz_to_rgb(xyz), spectrum_type)
    }

    /// Convert the SPD to XYZ cooefficients.
    fn to_xyz(&self) -> [Float; 3] {
        let (x, y, z) = (0..SPECTRAL_SAMPLES).fold((0.0, 0.0, 0.0), |(sx, sy, sz), i| {
            (
                sx + CIE_CURVES.x[i] * self.c[i],
                sy + CIE_CURVES.y[i] * self.c[i],
                sz + CIE_CURVES.z[i] * self.c[i],
            )
        });

        let scale = (SAMPLED_LAMBDA_END - SAMPLED_LAMBDA_START) as Float
            / (CIE_Y_INTEGRAL * SPECTRAL_SAMPLES as Float);

        [x * scale, y * scale, z * scale]
    }

    /// Returns the y-coefficient of XYZ colour.
    fn y(&self) -> Float {
        let yy = (0..SPECTRAL_SAMPLES).fold(0.0, |a, i| a + CIE_CURVES.y[i] * self.c[i]);
        yy * (SAMPLED_LAMBDA_END - SAMPLED_LAMBDA_START) as Float / (SPECTRAL_SAMPLES as Float)
    }

    /// Converts RGB values to a full SPD.
    ///
    /// * `rgb`           - RGB colour value.
    /// * `spectrum_type` - Indicates type of colour value. If `None`,
    ///                     defaults to `SpectrumType::Reflectance`.
    #[rustfmt::skip]
    fn from_rgb(rgb: &[Float; 3], spectrum_type: Option<SpectrumType>) -> Self {
        let spectrum_type = spectrum_type.map_or(SpectrumType::Reflectance, |s| s);
        let mut r = Self::default();

        // Convert spectrum to RGB.
        if rgb[0] <= rgb[1] && rgb[0] <= rgb[2] {
            // Compute reflectance `SampledSpectrum` with `rgb[0]` as minimum.
            r += rgb[0] * RGB_TO_SPECTRUM[spectrum_type][WHITE];
            if rgb[1] <= rgb[2] {
                r += (rgb[1] - rgb[0]) * RGB_TO_SPECTRUM[spectrum_type][CYAN];
                r += (rgb[2] - rgb[1]) * RGB_TO_SPECTRUM[spectrum_type][BLUE];
            } else {
                r += (rgb[2] - rgb[0]) * RGB_TO_SPECTRUM[spectrum_type][CYAN];
                r += (rgb[1] - rgb[2]) * RGB_TO_SPECTRUM[spectrum_type][GREEN];
            }
        } else if rgb[1] <= rgb[0] && rgb[1] <= rgb[2] {
            // Compute reflectance `SampledSpectrum` with 1rgb[1]1 as minimum.
            r += rgb[1] * RGB_TO_SPECTRUM[spectrum_type][WHITE];
            if rgb[0] <= rgb[2] {
                r += (rgb[0] - rgb[1]) * RGB_TO_SPECTRUM[spectrum_type][MAGENTA];
                r += (rgb[2] - rgb[0]) * RGB_TO_SPECTRUM[spectrum_type][BLUE];
            } else {
                r += (rgb[2] - rgb[1]) * RGB_TO_SPECTRUM[spectrum_type][MAGENTA];
                r += (rgb[0] - rgb[2]) * RGB_TO_SPECTRUM[spectrum_type][RED];
            }
        } else {
            // Compute reflectance `SampledSpectrum` with `rgb[2]` as minimum.
            r += rgb[2] * RGB_TO_SPECTRUM[spectrum_type][WHITE];
            if rgb[0] <= rgb[1] {
                r += (rgb[0] - rgb[2]) * RGB_TO_SPECTRUM[spectrum_type][YELLOW];
                r += (rgb[1] - rgb[0]) * RGB_TO_SPECTRUM[spectrum_type][GREEN];
            } else {
                r += (rgb[1] - rgb[2]) * RGB_TO_SPECTRUM[spectrum_type][YELLOW];
                r += (rgb[0] - rgb[1]) * RGB_TO_SPECTRUM[spectrum_type][RED];
            }
        }

        match spectrum_type {
            SpectrumType::Reflectance => r * 0.94,
            SpectrumType::Illuminant => r * 0.86445,
        }
    }

    /// Convert the SPD to RGB cooefficients.
    fn to_rgb(&self) -> [Float; 3] {
        xyz_to_rgb(&self.to_xyz())
    }

    /// Takes the square root of all sample values.
    fn sqrt(&self) -> Self {
        let c = self.c.map(|v| v.sqrt());
        Self { c }
    }

    /// Raises the sample values to a given power.
    ///
    /// * `p` - The power.
    fn pow(&self, p: Float) -> Self {
        let c = self.c.map(|v| v.powf(p));
        Self { c }
    }

    /// Converts to an `RGBSpectrum`.
    fn to_rgb_spectrum(&self) -> RGBSpectrum {
        RGBSpectrum::from(self.to_rgb())
    }

    /// Sets sample values `v` to `e^v`.
    ///
    /// * `other` - The other SPD.
    fn exp(&self) -> Self {
        let mut ret = self.clone();
        for s in ret.samples_mut().iter_mut() {
            *s = s.exp();
        }
        debug_assert!(!ret.has_nans());
        ret
    }
}

impl Add for SampledSpectrum {
    type Output = Self;

    /// Adds the corresponding sample values from another `SampledSpectrum`.
    ///
    /// * `other` - The other `SampledSpectrum`.
    fn add(self, other: Self) -> Self::Output {
        let mut ret = self;
        CoefficientSpectrum::add(&mut ret, &other);
        ret
    }
}

impl AddAssign for SampledSpectrum {
    /// Adds the corresponding sample values from another `SampledSpectrum`.
    ///
    /// * `other` - The other `SampledSpectrum`.
    fn add_assign(&mut self, other: Self) {
        CoefficientSpectrum::add(self, &other);
    }
}

impl Sub for SampledSpectrum {
    type Output = Self;

    /// Subtracts the corresponding sample values from another `SampledSpectrum`.
    ///
    /// * `other` - The other `SampledSpectrum`.
    fn sub(self, other: Self) -> Self::Output {
        let mut ret = self;
        CoefficientSpectrum::sub(&mut ret, &other);
        ret
    }
}

impl SubAssign for SampledSpectrum {
    /// Subtracts the corresponding sample values from another `SampledSpectrum`.
    ///
    /// * `other` - The other `SampledSpectrum`.
    fn sub_assign(&mut self, other: Self) {
        CoefficientSpectrum::sub(self, &other);
    }
}

impl Mul for SampledSpectrum {
    type Output = Self;

    /// Multiplies the corresponding sample values from another `SampledSpectrum`.
    ///
    /// * `other` - The other `SampledSpectrum`.
    fn mul(self, other: Self) -> Self::Output {
        let mut ret = self;
        CoefficientSpectrum::mul(&mut ret, &other);
        ret
    }
}

impl Mul<Float> for SampledSpectrum {
    type Output = Self;

    /// Scales the sample values with a constant factor.
    ///
    /// * `f` - Scaling factor.
    fn mul(self, f: Float) -> Self::Output {
        let mut ret = self;
        CoefficientSpectrum::scale(&mut ret, f);
        ret
    }
}

impl Mul<SampledSpectrum> for Float {
    type Output = SampledSpectrum;

    /// Scales the sample values of a `SampledSpectrum`.
    ///
    /// * `s` - Sample values.
    fn mul(self, s: SampledSpectrum) -> Self::Output {
        s * self
    }
}

impl Mul<&SampledSpectrum> for Float {
    type Output = SampledSpectrum;

    /// Scales the sample values of a `SampledSpectrum`.
    ///
    /// * `s` - Sample values.
    fn mul(self, s: &SampledSpectrum) -> Self::Output {
        *s * self
    }
}

impl MulAssign for SampledSpectrum {
    /// Multiplies the corresponding sample values from another `SampledSpectrum`.
    ///
    /// * `other` - The other `SampledSpectrum`.
    fn mul_assign(&mut self, other: Self) {
        CoefficientSpectrum::mul(self, &other);
    }
}

impl MulAssign<Float> for SampledSpectrum {
    /// Scales the sample values with a constant factor.
    ///
    /// * `f` - Scaling factor.
    fn mul_assign(&mut self, f: Float) {
        CoefficientSpectrum::scale(self, f);
    }
}

impl Div for SampledSpectrum {
    type Output = Self;

    /// Divides the corresponding sample values from another
    /// `SampledSpectrum`.
    ///
    /// * `other` - The other `SampledSpectrum`.
    fn div(self, other: Self) -> Self::Output {
        let mut ret = self;
        CoefficientSpectrum::div(&mut ret, &other);
        ret
    }
}

impl Div<Float> for SampledSpectrum {
    type Output = Self;

    /// Divides the sample values with given factor.
    ///
    /// * `f` - Scaling value.
    fn div(self, f: Float) -> Self::Output {
        let mut ret = self;
        CoefficientSpectrum::scale(&mut ret, 1.0 / f);
        ret
    }
}

impl DivAssign for SampledSpectrum {
    /// Divides the corresponding sample values from another `SampledSpectrum`.
    ///
    /// * `other` - The other `SampledSpectrum`.
    fn div_assign(&mut self, other: Self) {
        CoefficientSpectrum::div(self, &other);
    }
}

impl DivAssign<Float> for SampledSpectrum {
    /// Divides the sample values with given factor.
    ///
    /// * `f` - Scaling value.
    fn div_assign(&mut self, f: Float) {
        CoefficientSpectrum::scale(self, 1.0 / f);
    }
}

impl Neg for SampledSpectrum {
    type Output = Self;

    /// Scale the values by -1.
    fn neg(self) -> Self::Output {
        self * -1.0
    }
}

impl Index<usize> for SampledSpectrum {
    type Output = Float;

    /// Index the sample value.
    ///
    /// * `i` -  The index.
    fn index(&self, index: usize) -> &Self::Output {
        &self.c[index]
    }
}

impl IndexMut<usize> for SampledSpectrum {
    /// Index the sample to get a mutable sample value.
    ///
    /// * `i` - The index.
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        &mut self.c[i]
    }
}

impl Clamp<Float> for SampledSpectrum {
    /// Clamps the sample values.
    ///
    /// * `low`  - Low value.
    /// * `high` - High value.
    fn clamp(&self, low: Float, high: Float) -> Self {
        let c = self.c.map(|v| clamp(v, low, high));
        Self { c }
    }

    /// Clamps the values to [0.0, INFINITY].
    fn clamp_default(&self) -> Self {
        self.clamp(0.0, INFINITY)
    }
}

impl fmt::Display for SampledSpectrum {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[")?;
        for (i, v) in self.c.iter().enumerate() {
            write!(f, "{}", v)?;
            if i < SPECTRAL_SAMPLES - 1 {
                write!(f, ", ")?;
            }
        }
        write!(f, "]")
    }
}
