//! Utility functions.

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::pbrt::*;

/// Returns a direction (x, y, z) for spherical coordinates (θ, Ø).
///
/// * `sin_theta` - sin(θ).
/// * `cos_theta` - cos(θ).
/// * `phi`       - Ø.
#[inline]
pub fn spherical_direction(sin_theta: Float, cos_theta: Float, phi: Float) -> Vector3f {
    Vector3f::new(sin_theta * cos(phi), sin_theta * sin(phi), cos_theta)
}

/// Returns a direction (x, y, z) for spherical coordinates (θ, Ø) with respect
/// to a coordinate frame.
///
/// * `sin_theta` - sin(θ).
/// * `cos_theta` - cos(θ).
/// * `phi`       - Ø.
/// * `x`         - Basis vector representing x-axis.
/// * `y`         - Basis vector representing y-axis.
/// * `z`         - Basis vector representing z-axis.
#[inline]
pub fn spherical_direction_in_coord_frame(
    sin_theta: Float,
    cos_theta: Float,
    phi: Float,
    x: &Vector3f,
    y: &Vector3f,
    z: &Vector3f,
) -> Vector3f {
    sin_theta * cos(phi) * x + sin_theta * sin(phi) * y + cos_theta * z
}

/// Return the spherical angle θ for a given vector.
///
/// * `v` - The vector.
#[inline]
pub fn spherical_theta(v: &Vector3f) -> Float {
    clamp(v.z, -1.0, 1.0).acos()
}

/// Return the spherical angle Ø for a given vector.
///
/// * `v` - The vector.
#[inline]
pub fn spherical_phi(v: &Vector3f) -> Float {
    let p = atan2(v.y, v.x);
    if p < 0.0 {
        p + TWO_PI
    } else {
        p
    }
}
