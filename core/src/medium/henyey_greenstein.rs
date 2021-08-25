//! Henyey-Greenstein

use super::PhaseFunction;
use crate::geometry::*;
use crate::pbrt::*;

/// Henyey-Greenstein phase function.
struct HenyeyGreenstein {
    /// The asymmetry parameter. It is the average value of the product of the
    /// phase function being approximated and the cosine of the angle between two
    /// directions. Isotropic phase functions use g = 0.
    pub g: Float,
}
impl HenyeyGreenstein {
    /// Returns a new `HenyeyGreenstein`.
    ///
    /// * `g` - The asymmetry parameter.
    pub fn new(g: Float) -> Self {
        Self { g }
    }
}

impl PhaseFunction for HenyeyGreenstein {
    /// Returns the value of the phase function for the given pair of directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    fn p(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        phase_hg(wo.dot(wi), self.g)
    }

    /// Returns the phase function value and sampled incident direction given the
    /// outgoing direction and a sample value in [0, 1)^2.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - Sample value in [0, 1)^2.
    fn sample_p(&self, wo: &Vector3f, u: &Point2f) -> (Float, Vector3f) {
        // Compute $\cos \theta$ for Henyey--Greenstein sample
        let cos_theta = if abs(self.g) < 1e-3 {
            1.0 - 2.0 * u[0]
        } else {
            let sqr_term = (1.0 - self.g * self.g) / (1.0 + self.g - 2.0 * self.g * u[0]);
            -(1.0 + self.g * self.g - sqr_term * sqr_term) / (2.0 * self.g)
        };

        // Compute direction _wi_ for Henyey--Greenstein sample
        let sin_theta = max(0.0, 1.0 - cos_theta * cos_theta).sqrt();
        let phi = 2.0 * PI * u[1];

        let (v1, v2) = coordinate_system(&wo);
        let wi = spherical_direction_in_coord_frame(sin_theta, cos_theta, phi, &v1, &v2, wo);
        let phg = phase_hg(cos_theta, self.g);
        (phg, wi)
    }
}

/// Computes the Henyey-Greenstein phase function which can be used by other
/// phase function
///
/// * `cos_theta` - Angle between two direction vectors.
/// * `g`         - Asymmetry parametery.
#[inline]
pub fn phase_hg(cos_theta: Float, g: Float) -> Float {
    let denom = 1.0 + g * g + 2.0 * g * cos_theta;
    INV_FOUR_PI * (1.0 - g * g) / (denom * denom.sqrt())
}
