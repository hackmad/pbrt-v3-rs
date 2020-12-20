//! Triangle Filter

use super::{abs, max, Filter, FilterData, Float, Point2f, Vector2f};

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
    /// Returns value of the filter at a given point.
    ///
    /// * `p` - The position of the sample point relative to the center of the
    ///         filter. The point should be within the filter's extent.
    fn evaluate(&self, p: &Point2f) -> Float {
        max(0.0, self.data.radius.x - abs(p.x)) * max(0.0, self.data.radius.y - abs(p.y))
    }
}
