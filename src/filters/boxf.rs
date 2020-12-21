//! Box Filter

use super::{Filter, FilterData, Float, Point2f, Vector2f};

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
