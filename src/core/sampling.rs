//! Sampling functions

#![allow(dead_code)]
use super::geometry::*;
use super::pbrt::*;

/// Sample a point on a unit disk by mapping from a unit square to the unit
/// circle. The concentric mapping takes points in [-1, 1]^2 to unit disk by
/// uniformly mapping concentric squares to concentric circles.
///
/// * `u` - The random sample point.
pub fn concentric_sample_disk(u: &Point2f) -> Point2f {
    // Map uniform random numbers to [-1,1]^2.
    let u_offset = 2.0 * u - Vector2f::new(1.0, 1.0);

    // Handle degeneracy at the origin.
    if u_offset.x == 0.0 && u_offset.y == 0.0 {
        return Point2f::zero();
    }

    // Apply concentric mapping to point
    let (r, theta) = if abs(u_offset.x) > abs(u_offset.y) {
        (u_offset.x, PI_OVER_FOUR * (u_offset.y / u_offset.x))
    } else {
        (
            u_offset.y,
            PI_OVER_TWO - PI_OVER_FOUR * (u_offset.x / u_offset.y),
        )
    };

    r * Point2f::new(cos(theta), sin(theta))
}
