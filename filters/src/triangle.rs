//! Triangle Filter

use core::filter::*;
use core::geometry::*;
use core::paramset::*;
use core::pbrt::*;

/// Implements the triangle filter in which the weight falls off linearly from
/// the filter center over the square extent of the filter.
pub struct TriangleFilter {
    /// Filter data.
    pub data: FilterData,
}

impl TriangleFilter {
    /// Returns a new instance of `TriangleFilter`.
    ///
    /// * `radius` - Radius of the filter in x and y directions; beyond this
    ///              filter is 0.
    pub fn new(radius: Vector2f) -> Self {
        Self {
            data: FilterData::new(radius),
        }
    }
}

impl Filter for TriangleFilter {
    /// Return the filter parameters.
    fn get_data(&self) -> &FilterData {
        &self.data
    }

    /// Returns value of the filter at a given point.
    ///
    /// * `p` - The position of the sample point relative to the center of the
    ///         filter. The point should be within the filter's extent.
    fn evaluate(&self, p: &Point2f) -> Float {
        max(0.0, self.data.radius.x - abs(p.x)) * max(0.0, self.data.radius.y - abs(p.y))
    }
}

impl From<&ParamSet> for TriangleFilter {
    /// Create a `TriangleFilter` from `ParamSet`.
    ///
    /// * `params` - Parameter set.
    fn from(params: &ParamSet) -> Self {
        let xw = params.find_one_float("xwidth", 2.0);
        let yw = params.find_one_float("ywidth", 2.0);
        Self::new(Vector2f::new(xw, yw))
    }
}
