//! Mitchell Filter

#![allow(dead_code)]
use crate::core::filter::*;
use crate::core::geometry::*;
use crate::core::paramset::*;
use crate::core::pbrt::*;

/// Implements the Mitchell filter.
pub struct MitchellFilter {
    /// Filter data.
    pub data: FilterData,

    /// Parameter `B`.
    pub b: Float,

    /// Parameter `C`.
    pub c: Float,
}

impl MitchellFilter {
    /// Returns a new instance of `MitchellFilter`. Ideally the 2 parameters
    /// `B` and `C` should satisfy `B + 2C = 1`.
    ///
    /// * `radius` - Radius of the filter in x and y directions; beyond this
    ///              filter is 0.
    /// * `b`      - Parameter `B`.
    /// * `c`      - Parameter `C`.
    pub fn new(radius: Vector2f, b: Float, c: Float) -> Self {
        Self {
            data: FilterData::new(radius),
            b,
            c,
        }
    }

    /// Calculates the 1D filter function.
    ///
    /// * `x` - Distance from center of filter.
    fn mitchell_1d(&self, x: Float) -> Float {
        let x = abs(2.0 * x);

        if x > 1.0 {
            ((-self.b - 6.0 * self.c) * x * x * x
                + (6.0 * self.b + 30.0 * self.c) * x * x
                + (-12.0 * self.b - 48.0 * self.c) * x
                + (8.0 * self.c + 24.0 * self.c))
                * (1.0 / 6.0)
        } else {
            ((12.0 - 9.0 * self.b - 6.0 * self.c) * x * x * x
                + (-18.0 + 12.0 * self.b + 6.0 * self.c) * x * x
                + (6.0 - 2.0 * self.b))
                * (1.0 / 6.0)
        }
    }
}

impl Filter for MitchellFilter {
    /// Return the filter parameters.
    fn get_data(&self) -> &FilterData {
        &self.data
    }

    /// Returns value of the filter at a given point.
    ///
    /// * `p` - The position of the sample point relative to the center of the
    ///         filter. The point should be within the filter's extent.
    fn evaluate(&self, p: &Point2f) -> Float {
        self.mitchell_1d(p.x * self.data.inv_radius.x)
            * self.mitchell_1d(p.y * self.data.inv_radius.y)
    }
}

impl From<&ParamSet> for MitchellFilter {
    /// Create a `MitchellFilter` from `ParamSet`.
    ///
    /// * `params` - Parameter set.
    fn from(params: &ParamSet) -> Self {
        let xw = params.find_one_float("xwidth", 2.0);
        let yw = params.find_one_float("ywidth", 2.0);
        let b = params.find_one_float("B", 1.0 / 3.0);
        let c = params.find_one_float("C", 1.0 / 3.0);
        Self::new(Vector2f::new(xw, yw), b, c)
    }
}
