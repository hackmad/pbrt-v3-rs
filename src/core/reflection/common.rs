//! Common

#![allow(dead_code)]
use super::{abs, clamp, max, Dot, Float, Normal3f, Vector3f};

/// Returns the cosine of the angle θ measured from the given direction to the
/// z-axis.
///
/// * `w` - The direction vector.
#[inline]
pub fn cos_theta(w: &Vector3f) -> Float {
    w.z
}

/// Returns the square of the cosine of the angle θ measured from the given
/// direction to the z-axis.
///
/// * `w` - The direction vector.
#[inline]
pub fn cos_2_theta(w: &Vector3f) -> Float {
    w.z * w.z
}

/// Returns the absolute value of the cosine of the angle θ measured from the
/// given direction to the z-axis.
///
/// * `w` - The direction vector.
#[inline]
pub fn abs_cos_theta(w: &Vector3f) -> Float {
    abs(w.z)
}

/// Returns the square of the sine of the angle θ measured from the given
/// direction to the z-axis.
///
/// * `w` - The direction vector.
#[inline]
pub fn sin_2_theta(w: &Vector3f) -> Float {
    max(0.0, 1.0 - cos_2_theta(w))
}

/// Returns the sine of the angle θ measured from the given direction to the
/// z-axis.
///
/// * `w` - The direction vector.
#[inline]
pub fn sin_theta(w: &Vector3f) -> Float {
    sin_2_theta(w).sqrt()
}

/// Returns the tangent of the angle θ measured from the given direction to the
/// z-axis.
///
/// * `w` - The direction vector.
#[inline]
pub fn tan_theta(w: &Vector3f) -> Float {
    sin_theta(w) / cos_theta(w)
}

/// Returns the square of the tangent of the angle θ measured from the given
/// direction to the z-axis.
///
/// * `w` - The direction vector.
#[inline]
pub fn tan_2_theta(w: &Vector3f) -> Float {
    sin_2_theta(w) / cos_2_theta(w)
}

/// Returns the cosine of the angle Φ measured from the given direction to the
/// x-axis after projection to the xy plane.
///
/// * `w` - The direction vector.
#[inline]
pub fn cos_phi(w: &Vector3f) -> Float {
    let s = sin_theta(w);
    if s == 0.0 {
        1.0
    } else {
        clamp(w.x / s, -1.0, 1.0)
    }
}

/// Returns the square of the cosine of the angle Φ measured from the given
/// direction to the x-axis after projection to the xy plane.
///
/// * `w` - The direction vector.
#[inline]
pub fn cos_2_phi(w: &Vector3f) -> Float {
    let c = cos_phi(w);
    c * c
}

/// Returns the sine of the angle Φ measured from the given direction to the
/// x-axis after projection to the xy plane.
///
/// * `w` - The direction vector.
#[inline]
pub fn sin_phi(w: &Vector3f) -> Float {
    let s = sin_theta(w);
    if s == 0.0 {
        0.0
    } else {
        clamp(w.y / s, -1.0, 1.0)
    }
}

/// Returns the square of the sine of the angle Φ measured from the given
/// direction to the x-axis after projection to the xy plane.
///
/// * `w` - The direction vector.
#[inline]
pub fn sin_2_phi(w: &Vector3f) -> Float {
    let c = sin_phi(w);
    c * c
}

/// Returns the cosine of the angle ΔΦ between two vector's Φ values in the
/// shading coordinate system.
///
/// * `wa` - First direction vector.
/// * `wb` - Second direction vector.
pub fn cos_d_phi(wa: &Vector3f, wb: &Vector3f) -> Float {
    let waxy = wa.x * wa.x + wa.y * wa.y;
    let wbxy = wb.x * wb.x + wb.y * wb.y;
    if waxy == 0.0 || wbxy == 0.0 {
        1.0
    } else {
        clamp(
            (wa.x * wb.x + wa.y * wb.y) / (waxy * wbxy).sqrt(),
            -1.0,
            1.0,
        )
    }
}

/// Returns `true` if two vectors are in the same hemisphere.
///
/// * `w`  - First vector.
/// * `wp` - Second vector.
#[inline]
pub fn same_hemisphere(w: &Vector3f, wp: &Vector3f) -> bool {
    w.z * wp.z > 0.0
}

/// Computes the refracted direction, given incident direction `wi`, surface normal
/// `n` in the same hemisphere as `wi` and `eta`. If there is total internal
/// reflection, `None` is returned.
///
/// * `wi`  - Incident direction.
/// * `n`   - Surface normal.
/// * `eta` - Ratio of indices of refraction in the incident and transmitted media.
pub fn refract(wi: &Vector3f, n: &Normal3f, eta: Float) -> Option<Vector3f> {
    // Compute cos(theta_t) using Snell's law
    let cos_theta_i = n.dot(wi);
    let sin_2_theta_i = max(0.0, 1.0 - cos_theta_i * cos_theta_i);
    let sin_2_theta_t = eta * eta * sin_2_theta_i;

    // Handle total internal reflection for transmission.
    if sin_2_theta_t >= 1.0 {
        None
    } else {
        let cos_theta_t = (1.0 - sin_2_theta_t).sqrt();
        Some(eta * -(*wi) + (eta * cos_theta_i - cos_theta_t) * Vector3f::from(*n))
    }
}

/// Computes the reflection of a vector around a normal.
///
/// * `wo` - Vector to reflect.
/// * `n`  - Normal.
#[inline]
pub fn reflect(wo: &Vector3f, n: &Vector3f) -> Vector3f {
    -(*wo) + 2.0 * wo.dot(n) * n
}
