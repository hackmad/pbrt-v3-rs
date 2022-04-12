//! Common.

use crate::pbrt::*;
use crate::spectrum::RGBSpectrum;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
};

/// Determines if RGB value represents surface reflectance or illuminant.
#[derive(Copy, Clone)]
pub enum SpectrumType {
    /// Surface reflectance.
    Reflectance = 0,

    /// Surface illuminant.
    Illuminant = 1,
}
impl From<SpectrumType> for usize {
    fn from(stype: SpectrumType) -> Self {
        stype as usize
    }
}
pub const NUM_RGB_SPECTRUM_TYPES: usize = 2;

/// Specific colours in the RGB spectrum.
pub const WHITE: usize = 0;
pub const CYAN: usize = 1;
pub const MAGENTA: usize = 2;
pub const YELLOW: usize = 3;
pub const RED: usize = 4;
pub const GREEN: usize = 5;
pub const BLUE: usize = 6;
pub const NUM_RGB_SPECTRUM_COLOURS: usize = 7;

/// Stores a spectrum sample value at a given wavelenght.
#[derive(Copy, Clone, Default, Debug, PartialEq, PartialOrd)]
pub struct Sample {
    /// The wavelength.
    pub lambda: Float,

    /// The sample value.
    pub value: Float,
}
impl Sample {
    /// Create a new `Sample`.
    ///
    /// * `lambda` - The wavelength.
    /// * `value`  - The sample value.
    pub fn new(lambda: Float, value: Float) -> Self {
        Self { lambda, value }
    }

    /// Loads spectrum samples read from a data file.
    ///
    /// * `value` - A list of (wavelength, value) pairs. If list does not
    ///             contain even number of elements, the last one is ignored.
    pub fn list(values: &[Float]) -> Vec<Self> {
        let n = values.len();
        if n % 2 > 0 {
            warn!("Ignoring extra values in Sample::list().");
        }
        (0..n)
            .step_by(2)
            .map(|i| Sample::new(values[i], values[i + 1]))
            .collect()
    }
}

/// Interface and helper functions for SPDs.
pub trait CoefficientSpectrum:
    Sized
    + Add
    + AddAssign
    + Sub
    + SubAssign
    + Mul<Self>
    + MulAssign<Self>
    + MulAssign<Float>
    + Div<Self>
    + DivAssign<Self>
    + DivAssign<Float>
    + Neg
    + Index<usize>
    + IndexMut<usize>
    + Clamp<Float>
{
    /// Returns the stored samples.
    fn samples(&self) -> &[Float];

    /// Returns stored samples as mutable.
    fn samples_mut(&mut self) -> &mut [Float];

    /// Returns true if either coordinate is NaN.
    fn has_nans(&self) -> bool {
        for v in self.samples().iter() {
            if v.is_nan() {
                return true;
            }
        }
        false
    }

    /// Returns true if the values are zero everywhere.
    fn is_black(&self) -> bool {
        for v in self.samples().iter() {
            if *v != 0.0 {
                return false;
            }
        }
        true
    }

    /// Takes the square root of all sample values.
    fn sqrt(&self) -> Self;

    /// Raises the sample values to a given power.
    ///
    /// * `p` - The power.
    fn pow(&self, p: Float) -> Self;

    /// Returns the maximum sample value.
    fn max_component_value(&self) -> Float {
        let samples = self.samples();
        assert!(!samples.is_empty());
        samples[1..].iter().fold(samples[0], |m, v| max(m, *v))
    }

    /// Converts XYZ values to a full SPD.
    ///
    /// * `xyz`           - XYZ colour value.
    /// * `spectrum_type` - Indicates type of colour value. If `None`,
    ///                     defaults to `SpectrumType::Reflectance`.
    fn from_xyz(xyz: &[Float; 3], spectrum_type: Option<SpectrumType>) -> Self;

    /// Convert the SPD to XYZ cooefficients.
    fn to_xyz(&self) -> [Float; 3];

    /// Returns the y-coefficient of XYZ colour.
    fn y(&self) -> Float;

    /// Converts RGB values to a full SPD.
    ///
    /// * `rgb`           - RGB colour value.
    /// * `spectrum_type` - Indicates type of colour value. If `None`,
    ///                     defaults to `SpectrumType::Reflectance`.
    fn from_rgb(rgb: &[Float; 3], spectrum_type: Option<SpectrumType>) -> Self;

    /// Convert the SPD to RGB cooefficients.
    fn to_rgb(&self) -> [Float; 3];

    /// Converts to an `RGBSpectrum`.
    fn to_rgb_spectrum(&self) -> RGBSpectrum;

    /// Adds the sample values from another SPD.
    ///
    /// * `other` - The other SPD.
    fn add(&mut self, other: &Self) {
        debug_assert!(!self.has_nans());
        let samples = self.samples_mut();
        let other_samples = other.samples();
        let n = samples.len();
        assert!(n == other_samples.len());
        for i in 0..n {
            samples[i] += other_samples[i];
        }
    }

    /// Subtract the sample values from another SPD.
    ///
    /// * `other` - The other SPD.
    fn sub(&mut self, other: &Self) {
        debug_assert!(!self.has_nans());
        let samples = self.samples_mut();
        let other_samples = other.samples();
        let n = samples.len();
        assert!(n == other_samples.len());
        for i in 0..n {
            samples[i] -= other_samples[i];
        }
    }

    /// Multiplies the sample values from another SPD.
    ///
    /// * `other` - The other SPD.
    fn mul(&mut self, other: &Self) {
        debug_assert!(!self.has_nans());
        let samples = self.samples_mut();
        let other_samples = other.samples();
        let n = samples.len();
        assert!(n == other_samples.len());
        for i in 0..n {
            samples[i] *= other_samples[i];
        }
    }

    /// Divides the sample values from another SPD.
    ///
    /// * `other` - The other SPD.
    fn div(&mut self, other: &Self) {
        debug_assert!(!self.has_nans());
        let samples = self.samples_mut();
        let other_samples = other.samples();
        let n = samples.len();
        assert!(n == other_samples.len());
        for i in 0..n {
            samples[i] /= other_samples[i];
        }
    }

    /// Scales the sample values by a constant factor.
    ///
    /// * `f` - The factor.
    fn scale(&mut self, f: Float) {
        for s in self.samples_mut().iter_mut() {
            *s *= f;
        }
        debug_assert!(!self.has_nans());
    }

    /// Sets sample values `v` to `e^v`.
    ///
    /// * `other` - The other SPD.
    fn exp(&self) -> Self;
}

/// Determines if given vector containing wavelengths is sorted.
///
/// * `lambda` - Vector containing wavelengths.
pub fn are_spectrum_samples_sorted(samples: &[Sample]) -> bool {
    let n = samples.len();
    match n {
        0 => true,
        1 => true,
        _ => !samples
            .iter()
            .zip(samples[1..].iter())
            .any(|(s1, s2)| s1.lambda > s2.lambda),
    }
}

/// Sorts the sample values by wavelength.
///
/// * `lambda` - Vector containing wavelengths.
/// * `vals`   - Vector containing sample values corresponding to `lambda`.
pub fn sort_spectrum_samples(samples: &mut [Sample]) {
    samples.sort_by(|s1, s2| s1.lambda.partial_cmp(&s2.lambda).unwrap());
}

/// Returns the average value across segments that are partially or fully
/// within the range of wavelengths.
///
/// * `samples`      - The sample values.
/// * `lambda_start` - Starting wavelength.
/// * `lambda_end`   - Ending wavelength.
pub fn average_spectrum_samples(
    samples: &[Sample],
    lambda_start: Float,
    lambda_end: Float,
) -> Float {
    assert!(lambda_start < lambda_end);

    let n = samples.len();

    // Handle cases with out-of-bounds range or single sample only.
    if lambda_end <= samples[0].lambda {
        return samples[0].value;
    }

    if lambda_start >= samples[n - 1].lambda {
        return samples[n - 1].value;
    }

    if n == 1 {
        return samples[0].value;
    }

    let mut sum = 0.0;

    // Add contributions of constant segments before/after samples.
    if lambda_start < samples[0].lambda {
        sum += samples[0].value * (samples[0].lambda - lambda_start);
    }

    if lambda_end > samples[n - 1].lambda {
        sum += samples[n - 1].value * (lambda_end - samples[n - 1].lambda);
    }

    // Advance to first relevant wavelength segment
    let mut i = 0;
    while lambda_start > samples[i + 1].lambda {
        i += 1;
    }
    assert!(i + 1 < n);

    // Loop over wavelength sample segments and add contributions.
    let interp = |w: Float, i: usize| -> Float {
        lerp(
            (w - samples[i].lambda) / (samples[i + 1].lambda - samples[i].lambda),
            samples[i].value,
            samples[i + 1].value,
        )
    };

    loop {
        if i + 1 < n && lambda_end >= samples[i].lambda {
            break;
        }

        let seg_lambda_start = max(lambda_start, samples[i].lambda);
        let seg_lambda_end = min(lambda_end, samples[i + 1].lambda);

        sum += 0.5
            * (interp(seg_lambda_start, i) + interp(seg_lambda_end, i))
            * (seg_lambda_end - seg_lambda_start);

        i += 1;
    }
    sum / (lambda_end - lambda_start)
}

/// Returns the value of an SPD at a given wavelength by interpolating a possibly
/// irregularly sampled wavelengths/values by linearly interpolating between two
/// sample values that bracket the given wavelength.
///
/// * `samples` - The sample values.
/// * `l`       - Wavelength at which to interpolate SPD.
pub fn interpolate_spectrum_samples(samples: &[Sample], l: Float) -> Float {
    let n = samples.len();

    if l <= samples[0].lambda {
        return samples[0].value;
    }
    if l >= samples[n - 1].lambda {
        return samples[n - 1].value;
    }

    let offset = find_interval(n, |index| samples[index].lambda <= l);

    assert!(l >= samples[offset].lambda && l <= samples[offset + 1].lambda);

    let t = (l - samples[offset].lambda) / (samples[offset + 1].lambda - samples[offset].lambda);
    lerp(t, samples[offset].value, samples[offset + 1].value)
}

/// Converts the given XYZ coefficients to RGB coefficients using RGB spectra
/// defined for high-definition TVs.
///
/// * `xyz` - The XYZ coefficients.
#[rustfmt::skip]
pub fn xyz_to_rgb(xyz: &[Float; 3]) -> [Float; 3] {
    [
         3.240479 * xyz[0] - 1.537150 * xyz[1] - 0.498535 * xyz[2],
        -0.969256 * xyz[0] + 1.875991 * xyz[1] + 0.041556 * xyz[2],
         0.055648 * xyz[0] - 0.204043 * xyz[1] + 1.057311 * xyz[2],
    ]
}

/// Converts the given RGB coefficients to XYZ coefficients using RGB spectra
/// defined for high-definition TVs.
///
/// * `rgb` - The RGB coefficients.
#[rustfmt::skip]
pub fn rgb_to_xyz(rgb: &[Float; 3]) -> [Float; 3] {
    [
        0.412453 * rgb[0] + 0.357580 * rgb[1] + 0.180423 * rgb[2],
        0.212671 * rgb[0] + 0.715160 * rgb[1] + 0.072169 * rgb[2],
        0.019334 * rgb[0] + 0.119193 * rgb[1] + 0.950227 * rgb[2],
    ]
}

/// Returns the emitted radiance at a given temperature and wavelengths for a
/// blackbody (perfect emitter).
///
/// * `lambda` - Wavelengths in nanometers.
/// * `t`      - Temperature in Kelvin.
pub fn blackbody(lambda: &[Float], t: Float) -> Vec<Float> {
    if t <= 0.0 {
        return vec![0.0; lambda.len()];
    }

    const C: Float = 299792458.0;
    const H: Float = 6.62606957e-34;
    const KB: Float = 1.3806488e-23;

    let n = lambda.len();
    let mut le = vec![0.0; n];
    for i in 0..n {
        // Compute emitted radiance for blackbody at wavelength `lambda[i]`.
        let l = lambda[i] * 1e-9; // Convert nanometers -> meters.
        let lambda5 = (l * l) * (l * l) * l;
        le[i] = (2.0 * H * C * C) / (lambda5 * (((H * C) / (l * KB * t)).exp() - 1.0));
        assert!(!le[i].is_nan());
    }
    le
}

/// Returns the normalized emitted radiance at a given temperature and wavelengths
/// for a blackbody (perfect emitter) based on maximum blackbody radiance.
///
/// * `lambda` - Wavelengths in nanometers.
/// * `t`      - Temperature in Kelvin.
pub fn blackbody_normalized(lambda: &[Float], t: Float) -> Vec<Float> {
    let mut le = blackbody(lambda, t);

    // Normalize `Le` values based on maximum blackbody radiance.
    let lambda_max = 2.8977721e-3 / t * 1e9; // Convert to meters -> nanometers.
    let max_l = blackbody(&[lambda_max], t);

    for v in &mut le {
        *v /= max_l[0];
    }
    le
}
