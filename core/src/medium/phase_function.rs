//! Phase Function

use super::HenyeyGreenstein;
use crate::geometry::*;
use crate::pbrt::*;

/// Models scattering properties in volumetric media.
pub enum PhaseFunction<'arena> {
    HenyeyGreenstein(&'arena mut HenyeyGreenstein),
}

impl<'arena> PhaseFunction<'arena> {
    /// Returns the value of the phase function for the given pair of directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    pub fn p(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        match self {
            PhaseFunction::HenyeyGreenstein(f) => f.p(wo, wi),
        }
    }

    /// Returns the phase function value and sampled incident direction given the
    /// outgoing direction and a sample value in [0, 1)^2.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - Sample value in [0, 1)^2.
    pub fn sample_p(&self, wo: &Vector3f, u: &Point2f) -> (Float, Vector3f) {
        match self {
            PhaseFunction::HenyeyGreenstein(f) => f.sample_p(wo, u),
        }
    }
}
