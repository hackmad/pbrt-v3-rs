//! 1D Distribution.

use crate::pbrt::*;

/// Represents a piecewise-constant 1D function’s PDF and CDF and provides methods to perform this sampling efficiently.
#[derive(Clone)]
pub struct Distribution1D {
    /// Piecewise-constant function.
    pub func: Vec<Float>,

    /// CDF for `func`.
    pub cdf: Vec<Float>,

    /// Integral of `func`.
    pub func_int: Float,
}

impl Distribution1D {
    /// Returns a new `Distribution1D` for given piecewise-constant function.
    ///
    /// - `f` - Piecewise-constant 1D function.
    pub fn new(f: Vec<Float>) -> Self {
        let n = f.len();

        // Compute integral of step function at `x_i`
        let mut cdf: Vec<Float> = Vec::with_capacity(n + 1);
        cdf.push(0.0);
        for i in 1..n + 1 {
            cdf.push(cdf[i - 1] + f[i - 1] / n as Float);
        }

        // Transform step function integral into CDF.
        let func_int = cdf[n];
        if func_int == 0.0 {
            for (i, v) in cdf.iter_mut().enumerate().skip(1).take(n) {
                *v = i as Float / n as Float;
            }
        } else {
            for v in cdf.iter_mut().skip(1).take(n) {
                *v /= func_int;
            }
        }

        Self { func: f, cdf, func_int }
    }

    /// Returns the number of sample points for the piecewise-constant function.
    pub fn count(&self) -> usize {
        self.func.len()
    }

    /// Return a sample in [0, 1), PDF and offset from the distribution given a random sample.
    ///
    /// - `u` - The random sample.
    pub fn sample_continuous(&self, u: Float) -> (Float, Float, usize) {
        // Find surrounding CDF segments and `offset`.
        let offset = find_interval(self.cdf.len(), |index| self.cdf[index] <= u);

        // Compute offset along CDF segment.
        let mut du = u - self.cdf[offset];
        if self.cdf[offset + 1] - self.cdf[offset] > 0.0 {
            assert!(self.cdf[offset + 1] > self.cdf[offset]);
            du /= self.cdf[offset + 1] - self.cdf[offset];
        }
        debug_assert!(!du.is_nan());

        // Compute PDF for sampled offset.
        let pdf = if self.func_int > 0.0 {
            self.func[offset] / self.func_int
        } else {
            0.0
        };

        // Return `x` in [0,1) corresponding to sample, PDF and offset.
        ((offset as Float + du) / self.count() as Float, pdf, offset)
    }

    /// Return a sample from the discrete distribution given a random sample.
    ///
    /// - `u` - The random sample.
    pub fn sample_discrete(&self, u: Float) -> (usize, Float, Float) {
        // Find surrounding CDF segments and `offset`.
        let offset = find_interval(self.cdf.len(), |index| self.cdf[index] <= u);
        let pdf = if self.func_int > 0.0 {
            self.func[offset] / (self.func_int * self.count() as Float)
        } else {
            0.0
        };
        let u_remapped = (u - self.cdf[offset]) / (self.cdf[offset + 1] - self.cdf[offset]);

        //assert!(u_remapped >= 0.0 && u_remapped <= 1.0);
        assert!((0.0..=1.0).contains(&u_remapped));

        (offset, pdf, u_remapped)
    }

    /// Return the PDF for sampling a given value from the discrete PDF.
    ///
    /// * `index` - Sample index.
    pub fn discrete_pdf(&self, index: usize) -> Float {
        assert!(index < self.count());
        self.func[index] / (self.func_int * self.count() as Float)
    }
}
