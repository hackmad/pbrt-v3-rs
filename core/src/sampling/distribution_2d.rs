//! 2D Distribution.

use crate::geometry::*;
use crate::pbrt::*;
use crate::sampling::Distribution1D;

/// Represents a piecewise-constant 2D function’s PDF and CDF and provides methods to perform this sampling efficiently.
#[derive(Clone)]
pub struct Distribution2D {
    /// 1D conditional sampling density `p[ũ|ṽ]` for each `nv`.
    p_conditional_v: Vec<Distribution1D>,

    /// Marginal sampling density p[ṽ].
    p_marginal: Distribution1D,
}

impl Distribution2D {
    /// Returns a new `Distribution2D` for given piecewise-constant function.
    ///
    /// - `func` - Piecewise-constant 2D function.
    pub fn new(func: Vec<Vec<Float>>) -> Self {
        let p_conditional_v: Vec<Distribution1D> = func.iter().map(|f| Distribution1D::new(f.clone())).collect();
        let marginal_func: Vec<Float> = p_conditional_v.iter().map(|pcv| pcv.func_int).collect();
        let p_marginal = Distribution1D::new(marginal_func);
        Self {
            p_conditional_v,
            p_marginal,
        }
    }

    /// Return a sample point and PDF from the distribution given a random sample.
    ///
    /// - `u` - The random sample.
    pub fn sample_continuous(&self, u: &Point2f) -> (Point2f, Float) {
        // Draw a sample from the p[ṽ] marginal distribution in order to find the
        // ṽ coordinate.
        let (d1, pdf1, v) = self.p_marginal.sample_continuous(u[1]);

        // Use ṽ to find the precomputed conditional distribution to use for
        // sampling ũ.
        let (d0, pdf0, _) = self.p_conditional_v[v].sample_continuous(u[0]);

        let pdf = pdf0 * pdf1;

        (Point2f::new(d0, d1), pdf)
    }

    /// Return the PDF value for a given sample value.
    ///
    /// * `p` - Sample value.
    pub fn pdf(&self, p: &Point2f) -> Float {
        // Compute the product of the conditional and marginal PDFs for sampling
        // it from the distribution.
        let iu = clamp(
            (p[0] * self.p_conditional_v[0].count() as Float) as usize,
            0_usize,
            self.p_conditional_v[0].count() - 1,
        );
        let iv = clamp(
            (p[1] * self.p_marginal.count() as Float) as usize,
            0_usize,
            self.p_marginal.count() - 1,
        );
        self.p_conditional_v[iv].func[iu] / self.p_marginal.func_int
    }
}
