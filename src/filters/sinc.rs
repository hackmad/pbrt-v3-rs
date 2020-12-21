//! LanczosSinc Filter

use super::{abs, sin, Filter, FilterData, Float, Point2f, Vector2f, PI};

/// Implements a windowed sinc filter.
pub struct LanczosSincFilter {
    /// Filter data.
    pub data: FilterData,

    /// Number of cycles the sinc function passes through before it is clamped
    /// to 0.
    pub tau: Float,
}

impl LanczosSincFilter {
    /// Returns a new instance of `LanczosSincFilter`.
    ///
    /// * `radius` - Radius of the filter in x and y directions; beyond this
    ///              filter is 0.
    /// * `tau`    - Number of cycles the sinc function passes through before
    ///              it is clamped to 0.
    pub fn new(radius: Vector2f, tau: Float) -> Self {
        Self {
            data: FilterData::new(radius),
            tau,
        }
    }

    /// Calculates the LanczosSinc filter function for a given distance.
    ///
    /// * `x`      - Distance from center of filter.
    /// * `radius` - Radius of filter beyond which it is 0.
    fn windowed_sinc(&self, x: Float, radius: Float) -> Float {
        let x = abs(x);
        if x > radius {
            0.0
        } else {
            let lanczos = sinc(x / self.tau);
            sinc(x) * lanczos
        }
    }
}

impl Filter for LanczosSincFilter {
    /// Return the filter parameters.
    fn get_data(&self) -> &FilterData {
        &self.data
    }

    /// Returns value of the filter at a given point.
    ///
    /// * `p` - The position of the sample point relative to the center of the
    ///         filter. The point should be within the filter's extent.
    fn evaluate(&self, p: &Point2f) -> Float {
        self.windowed_sinc(p.x, self.data.radius.x) * self.windowed_sinc(p.y, self.data.radius.y)
    }
}

/// Evaluates the sinc function.
///
/// `x` - Point to evaluate sinc function at.
fn sinc(x: Float) -> Float {
    let x = abs(x);
    if x < 1e-5 {
        1.0
    } else {
        sin(PI * x) / (PI * x)
    }
}
