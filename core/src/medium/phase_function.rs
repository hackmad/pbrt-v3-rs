//! Phase Function

use super::HenyeyGreenstein;
use crate::geometry::*;
use crate::pbrt::*;
use std::fmt;

/// Models scattering properties in volumetric media.
#[derive(Clone)]
pub enum PhaseFunction {
    HenyeyGreenstein(HenyeyGreenstein),
}

impl PhaseFunction {
    /// Returns the value of the phase function for the given pair of directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    pub fn p(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        match self {
            PhaseFunction::HenyeyGreenstein(f) => f.p(wo, wi),
        }
    }

    /// Returns the phase function value and sampled incident direction given the outgoing direction and a sample value
    /// in [0, 1)^2.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - Sample value in [0, 1)^2.
    pub fn sample_p(&self, wo: &Vector3f, u: &Point2f) -> (Float, Vector3f) {
        match self {
            PhaseFunction::HenyeyGreenstein(f) => f.sample_p(wo, u),
        }
    }
}

impl fmt::Display for PhaseFunction {
    /// Formats the value using the given formatter.
    ///
    /// * `f` - Formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[PhaseFunction ")?;
        match self {
            Self::HenyeyGreenstein(h) => write!(f, "{}", h)?,
        }
        write!(f, "]")
    }
}
