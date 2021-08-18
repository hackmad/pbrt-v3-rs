//! Interval of real numbers

#![allow(dead_code)]
use crate::core::pbrt::*;
use std::fmt;
use std::mem::swap;
use std::ops::{Add, Mul, Sub};

/// Interval of real numbers used to give conservative bounds for functions.
#[derive(Copy, Clone, Default)]
pub struct Interval {
    // Low value.
    pub low: Float,

    // High value.
    pub high: Float,
}

impl Interval {
    /// Create an interval with given bounds. The interval will sort the input.
    ///
    /// * `v0` - A real number.
    /// * `v1` - A real number.
    pub fn new(v0: Float, v1: Float) -> Self {
        Self {
            low: min(v0, v1),
            high: max(v0, v1),
        }
    }

    /// Return the interval used for sine function assuming this interval fits
    /// inside [0, 2π].
    pub fn sin(&self) -> Interval {
        assert!(
            self.low >= 0.0,
            "interval low < 0 not allowed for sine function"
        );
        assert!(
            self.high <= 2.0001 * PI,
            "interval high > 2π not allowed for sine function"
        );

        let mut sin_low = self.low.sin();
        let mut sin_high = self.high.sin();

        if sin_low > sin_high {
            swap(&mut sin_low, &mut sin_high);
        }

        if self.low < PI_OVER_TWO && self.high > PI_OVER_TWO {
            sin_high = 1.0;
        }

        if self.low < 3.0 * PI_OVER_TWO && self.high > 3.0 * PI_OVER_TWO {
            sin_low = -1.0;
        }

        Interval::new(sin_low, sin_high)
    }

    /// Return the interval used for cosine function assuming this interval fits
    /// inside [0, 2π].
    pub fn cos(&self) -> Interval {
        assert!(
            self.low >= 0.0,
            "interval low < 0 not allowed for cosine function"
        );
        assert!(
            self.high <= 2.0001 * PI,
            "interval high > 2π not allowed for cosine function"
        );

        let mut cos_low = self.low.cos();
        let mut cos_high = self.high.cos();

        if cos_low > cos_high {
            swap(&mut cos_low, &mut cos_high);
        }

        if self.low < PI && self.high > PI {
            cos_low = -1.0;
        }

        Interval::new(cos_low, cos_high)
    }

    /// Find the values of any zero crossings of the equation:
    ///
    /// d/dt(a_M,p(t)) = c1 + (c2 + c3 * t) * cos(2θ * t) + (c4 + c5 * t) * sin(2θ * t)
    ///
    /// over the interval.
    ///
    /// * `c1`         - Coefficient c1.
    /// * `c2`         - Coefficient c2.
    /// * `c3`         - Coefficient c3.
    /// * `c4`         - Coefficient c4.
    /// * `c5`         - Coefficient c5.
    /// * `theta`      - Angle θ.
    /// * `zeros`      - Slice to return zeros in.
    /// * `zero_count` - Number of zeros found.
    /// * `depth`      - Depth of search (first call should be `zeros.len()`).
    pub fn find_zeros(
        &self,
        c1: Float,
        c2: Float,
        c3: Float,
        c4: Float,
        c5: Float,
        theta: Float,
        zeros: &mut [Float; 8],
        zero_count: &mut usize,
        depth: usize,
    ) {
        // Evaluate motion derivative in interval form, return if no zeros
        let range = Interval::from(c1)
            + (Interval::from(c2) + Interval::from(c3) * *self)
                * (Interval::from(2.0 * theta) * *self).cos()
            + (Interval::from(c4) + Interval::from(c5) * *self)
                * (Interval::from(2.0 * theta) * *self).sin();

        if range.low > 0.0 || range.high < 0.0 || range.low == range.high {
            return;
        }

        if depth > 0 {
            // Split self and check both resulting intervals
            let mid = (self.low + self.high) * 0.5;

            Interval::new(self.low, mid).find_zeros(
                c1,
                c2,
                c3,
                c4,
                c5,
                theta,
                zeros,
                zero_count,
                depth - 1,
            );

            Interval::new(mid, self.high).find_zeros(
                c1,
                c2,
                c3,
                c4,
                c5,
                theta,
                zeros,
                zero_count,
                depth - 1,
            );
        } else {
            // Use Newton's method to refine zero
            let mut t_newton = (self.low + self.high) * 0.5;
            for _ in 0..4 {
                let f_newton = c1
                    + (c2 + c3 * t_newton) * (2.0 * theta * t_newton).cos()
                    + (c4 + c5 * t_newton) * (2.0 * theta * t_newton).sin();

                let f_prime_newton = (c3 + 2.0 * (c4 + c5 * t_newton) * theta)
                    * (2.0 * t_newton * theta).cos()
                    + (c5 - 2.0 * (c2 + c3 * t_newton) * theta) * (2.0 * t_newton * theta).sin();

                if f_newton == 0.0 || f_prime_newton == 0.0 {
                    break;
                }

                t_newton -= f_newton / f_prime_newton;
            }
            if t_newton >= self.low - 1e-3 && t_newton < self.high + 1e-3 {
                zeros[*zero_count] = t_newton;
                *zero_count += 1;
            }
        }
    }
}

impl From<Float> for Interval {
    /// Create an interval on single point.
    ///
    /// * `v` - A real number.
    fn from(v: Float) -> Self {
        Self { low: v, high: v }
    }
}

impl Add for Interval {
    type Output = Self;

    /// Returns the conservative bounds for addition.
    ///
    /// * `i` -  The interval to add.
    fn add(self, i: Self) -> Self::Output {
        Interval::new(self.low + i.low, self.high + i.high)
    }
}

impl Sub for Interval {
    type Output = Self;

    /// Returns the conservative bounds for subtraction.
    ///
    /// * `i` -  The interval to subtract.
    fn sub(self, i: Self) -> Self::Output {
        Interval::new(self.low - i.low, self.high - i.high)
    }
}

impl Mul for Interval {
    type Output = Self;

    /// Returns the conservative bounds for multiplication.
    ///
    /// * `i` -  The interval to multiply.
    fn mul(self, i: Self) -> Self::Output {
        // The sides of each interval which determine the minimum and maximum
        // values of the result interval depend on the signs of the respective
        // values. Multiplying the various possibilities and taking the
        // overall minimum and maximum is easier than working through which ones
        // to use and multiplying these.
        Self::new(
            min(
                min(self.low * i.low, self.high * i.low),
                min(self.low * i.high, self.high * i.high),
            ),
            max(
                max(self.low * i.low, self.high * i.low),
                max(self.low * i.high, self.high * i.high),
            ),
        )
    }
}

impl fmt::Debug for Interval {
    /// Formats the value using the given formatter.
    ///
    /// * `f` - Formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Interval")
            .field("low", &self.low)
            .field("high", &self.high)
            .finish()
    }
}

impl fmt::Display for Interval {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}, {}]", self.low, self.high)
    }
}
