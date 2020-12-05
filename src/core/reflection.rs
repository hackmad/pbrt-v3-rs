//! Reflection and surface scattering models.

use super::geometry::{Dot, Normal3f, Vector3f};
use super::pbrt::{max, Float};
use std::sync::Arc;

/// BSDF model.
pub struct BSDF {}

/// Atomic reference counted `BSDF`.
pub type ArcBSDF = Arc<BSDF>;

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
