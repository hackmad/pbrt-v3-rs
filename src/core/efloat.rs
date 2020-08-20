//! EFloat

#![allow(dead_code)]
use super::pbrt::*;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Same as regular float but uses operator overloading to provide
/// arithmetic operations while computing error bounds.
#[derive(Copy, Clone, Debug, Default)]
pub struct EFloat {
    /// The floating point value.
    v: f32,

    /// The lower bound on `v`.
    low: f32,

    /// The upper bound on `v`.
    high: f32,

    /// 64-bit precision value corresponding to `v`.
    v_precise: f64,
}

/// Construct a new float with error bounds.
///
/// * `v`   - The 32-bit floating point value.
/// * `err` - The error (default to 0.0).
pub fn efloat(v: f32, err: f32) -> EFloat {
    let mut r = EFloat::default();

    if err == 0.0 {
        r.low = v;
        r.high = v;
    } else {
        // Compute conservative bounds by rounding the endpoints away from the
        // middle. Note that this will be over-conservative in cases where v-err
        // or v+err are exactly representable in floating-point, but it's
        // probably not worth the trouble of checking this case.
        r.low = next_float_down(r.v - err);
        r.high = next_float_up(r.v + err);
    };

    r.v_precise = v as f64;

    r.check();
    r
}

/// Construct a new float with error bounds and 64-bit value.
///
/// * `v`   - The 32-bit floating point value.
/// * `ld`  - The 64-bit floating point value.
/// * `err` - The error (default to 0.0).
pub fn efloat_precise(v: f32, ld: f64, err: f32) -> EFloat {
    let mut r = efloat(v, err);
    r.v_precise = ld;
    r.check();
    r
}

/// Solve the quadratic equation a * x ^ 2  + b * x + c = 0.
///
/// * `a` - Coefficient of x ^ 2 term.
/// * `b` - Coefficient of x term.
/// * `c` - Coefficient of constant term.
pub fn quadratic(a: EFloat, b: EFloat, c: EFloat) -> Option<(EFloat, EFloat)> {
    // Find quadratic discriminant
    let discrim: f64 = b.v as f64 * b.v as f64 - 4.0f64 * a.v as f64 * c.v as f64;
    if discrim < 0.0 {
        None
    } else {
        let root_discrim = discrim.sqrt() as f32;
        let ef_root_discrim = efloat(root_discrim, MACHINE_EPSILON * root_discrim);

        // Compute quadratic _t_ values
        let q = if b.v < 0.0 {
            -0.5f32 * (b - ef_root_discrim)
        } else {
            -0.5f32 * (b + ef_root_discrim)
        };

        let t0 = q / a;
        let t1 = c / q;

        if t0.v > t1.v {
            Some((t1, t0))
        } else {
            Some((t0, t1))
        }
    }
}

impl EFloat {
    /// Asserts low < high, low < v_precise < high for finite non-NAN values.
    fn check(&self) {
        if self.low.is_finite()
            && !self.low.is_nan()
            && self.high.is_finite()
            && !self.high.is_nan()
        {
            assert_eq!(self.low < self.high, true);
        }

        if self.v.is_finite() && !self.v.is_nan() {
            assert_eq!((self.low as f64) < self.v_precise, true);
            assert_eq!(self.v_precise < (self.high as f64), true);
        }
    }

    /// Returns the lower bound on the original value.
    pub fn lower_bound(&self) -> f32 {
        self.low
    }

    /// Returns the upper bound on the original value.
    pub fn upper_bound(&self) -> f32 {
        self.high
    }

    /// Returns the 64-bit precision value.
    pub fn precise_value(&self) -> f64 {
        self.v_precise
    }

    /// Returns the relative error in the 64-bit precision value.
    pub fn relative_error(&self) -> f32 {
        (((self.v_precise - (self.v as f64)) / self.v_precise).abs()) as f32
    }

    /// Returns the absolute error.
    pub fn get_absolute_error(&self) -> f32 {
        next_float_up(max((self.high - self.v).abs(), (self.v - self.low).abs()))
    }

    /// Returns the square root.
    pub fn sqrt(&self, ef: EFloat) -> EFloat {
        let mut r = EFloat::default();

        r.v = ef.v.sqrt();
        r.v_precise = ef.v_precise.sqrt();

        r.low = next_float_down(ef.low.sqrt());
        r.high = next_float_up(ef.high.sqrt());

        r.check();
        r
    }

    /// Returns the absolute value.
    pub fn abs(&self) -> EFloat {
        if self.low >= 0.0 {
            // The entire interval is greater than zero, so we're all set.
            *self
        } else if self.high <= 0.0 {
            // The entire interval is less than zero.
            let mut r = EFloat::default();

            r.v = -self.v;
            r.v_precise = -self.v_precise;

            r.low = -self.high;
            r.high = -self.low;

            r.check();
            r
        } else {
            // The interval straddles zero.
            let mut r = EFloat::default();

            r.v = self.v.abs();
            r.v_precise = self.v_precise.abs();

            r.low = 0.0;
            r.high = max(-self.low, self.high);

            r.check();
            r
        }
    }
}

impl From<f32> for EFloat {
    /// Converts a 32-bit floating point value to EFloat.
    ///
    /// * `v` - The 32-bit floating point value to convert.
    fn from(v: f32) -> Self {
        efloat(v, 0.0)
    }
}

impl From<EFloat> for f32 {
    /// Converts EFloat to 32-bit floating point value.
    ///
    /// * `ef` - The EFloat value.
    fn from(ef: EFloat) -> f32 {
        ef.v
    }
}

impl From<EFloat> for f64 {
    /// Converts EFloat to 64-bit floating point value.
    ///
    /// * `ef` - The EFloat value.
    fn from(ef: EFloat) -> f64 {
        ef.v as f64
    }
}

impl PartialEq for EFloat {
    /// Implement partial equality based on the stored value `v`.
    ///
    /// * `other` - The other EFloat.
    fn eq(&self, other: &Self) -> bool {
        self.v == other.v
    }
}

impl Add for EFloat {
    type Output = Self;

    /// Add an EFloat.
    ///
    /// * `ef` - The value to add.
    fn add(self, ef: EFloat) -> Self::Output {
        let mut r = EFloat::default();

        r.v = self.v + ef.v;
        r.v_precise = self.v_precise + ef.v_precise;

        // Interval arithemetic addition, with the result rounded away from
        // the value r.v in order to be conservative.
        r.low = next_float_down(self.low + ef.low);
        r.high = next_float_up(self.high + ef.high);

        r.check();
        r
    }
}

impl Add<f32> for EFloat {
    type Output = Self;

    /// Add an f32.
    ///
    /// * `v` - The value to add.
    fn add(self, v: f32) -> Self::Output {
        self + EFloat::from(v)
    }
}

impl Add<EFloat> for f32 {
    type Output = EFloat;

    /// Add an EFloat.
    ///
    /// * `ef` - The value to add.
    fn add(self, ef: EFloat) -> Self::Output {
        EFloat::from(self) + ef
    }
}

impl Sub for EFloat {
    type Output = Self;

    /// Subtract an EFloat.
    ///
    /// * `ef` - The value to subtract.
    fn sub(self, ef: EFloat) -> Self::Output {
        let mut r = EFloat::default();

        r.v = self.v - ef.v;
        r.v_precise = self.v_precise - ef.v_precise;

        r.low = next_float_down(self.low - ef.high);
        r.high = next_float_up(self.high - ef.low);

        r.check();
        r
    }
}

impl Sub<f32> for EFloat {
    type Output = Self;

    /// Subtract an f32.
    ///
    /// * `v` - The value to subtract.
    fn sub(self, v: f32) -> Self::Output {
        self - EFloat::from(v)
    }
}

impl Sub<EFloat> for f32 {
    type Output = EFloat;

    /// Subtract an EFloat.
    ///
    /// * `ef` - The value to subtract
    fn sub(self, ef: EFloat) -> Self::Output {
        EFloat::from(self) - ef
    }
}

impl Mul for EFloat {
    type Output = Self;

    /// Multiply an EFloat.
    ///
    /// * `ef` - The value to multiply.
    fn mul(self, ef: EFloat) -> Self::Output {
        let mut r = EFloat::default();

        r.v = self.v * ef.v;
        r.v_precise = self.v_precise * ef.v_precise;

        let prod = [
            self.low * ef.low,
            self.high * ef.low,
            self.low * ef.high,
            self.high * ef.high,
        ];

        r.low = next_float_down(min(min(prod[0], prod[1]), min(prod[2], prod[3])));
        r.high = next_float_up(max(max(prod[0], prod[1]), max(prod[2], prod[3])));

        r.check();
        r
    }
}

impl Mul<f32> for EFloat {
    type Output = Self;

    /// Multiply an f32.
    ///
    /// * `v` - The value to multiply.
    fn mul(self, v: f32) -> Self::Output {
        self * EFloat::from(v)
    }
}

impl Mul<EFloat> for f32 {
    type Output = EFloat;

    /// Multiply an EFloat.
    ///
    /// * `ef` - The value to multiply.
    fn mul(self, ef: EFloat) -> Self::Output {
        EFloat::from(self) * ef
    }
}

impl Div for EFloat {
    type Output = Self;

    /// Divide by an EFloat.
    ///
    /// * `ef` - The value to divide by.
    fn div(self, ef: EFloat) -> Self::Output {
        let mut r = EFloat::default();

        r.v = self.v / ef.v;
        r.v_precise = self.v_precise / ef.v_precise;

        if ef.low < 0.0 && ef.high > 0.0 {
            // The interval we're dividing by straddles zero, so just
            // return an interval of everything.
            r.low = f32::NEG_INFINITY;
            r.high = f32::INFINITY;
        } else {
            let div = [
                self.low / ef.low,
                self.high / ef.low,
                self.low / ef.high,
                self.high / ef.high,
            ];
            r.low = next_float_down(min(min(div[0], div[1]), min(div[2], div[3])));
            r.high = next_float_up(max(max(div[0], div[1]), max(div[2], div[3])));
        }

        r.check();
        r
    }
}

impl Div<f32> for EFloat {
    type Output = Self;

    /// Divide by f32.
    ///
    /// * `v` - The value to divide by.
    fn div(self, v: f32) -> Self::Output {
        self / EFloat::from(v)
    }
}

impl Div<EFloat> for f32 {
    type Output = EFloat;

    /// Divide by EFloat.
    ///
    /// * `ef` - The value to divide by.
    fn div(self, ef: EFloat) -> Self::Output {
        EFloat::from(self) / ef
    }
}

impl Neg for EFloat {
    type Output = Self;

    /// Return the negative value.
    fn neg(self) -> Self::Output {
        let mut r = EFloat::default();

        r.v = -self.v;
        r.v_precise = -self.v_precise;

        r.low = -self.high;
        r.high = -self.low;

        r.check();
        r
    }
}
