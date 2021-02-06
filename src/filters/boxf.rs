//! Box Filter

#![allow(dead_code)]
use crate::core::filter::*;
use crate::core::geometry::*;
use crate::core::paramset::*;
use crate::core::pbrt::*;

/// Implements the box filter which equally weights all samples within a square
/// region of the image.
pub struct BoxFilter {
    /// Filter data.
    pub data: FilterData,
}

impl BoxFilter {
    /// Returns a new instance of `BoxFilter`.
    ///
    /// * `radius` - Radius of the filter in x and y directions; beyond this
    ///              filter is 0.
    pub fn new(radius: Vector2f) -> Self {
        Self {
            data: FilterData::new(radius),
        }
    }
}

impl Filter for BoxFilter {
    /// Return the filter parameters.
    fn get_data(&self) -> &FilterData {
        &self.data
    }

    /// Returns value of the filter at a given point.
    ///
    /// * `p` - The position of the sample point relative to the center of the
    ///         filter. The point should be within the filter's extent.
    fn evaluate(&self, _p: &Point2f) -> Float {
        1.0
    }
}

impl From<&ParamSet> for BoxFilter {
    /// Create a `BoxFilter` from `ParamSet`.
    ///
    /// * `params` - Parameter set.
    fn from(params: &ParamSet) -> Self {
        let xw = params.find_one_float("xwidth", 0.5);
        let yw = params.find_one_float("ywidth", 0.5);
        Self::new(Vector2f::new(xw, yw))
    }
}
