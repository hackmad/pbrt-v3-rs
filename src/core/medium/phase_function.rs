//! Phase Function

use crate::core::geometry::*;
use crate::core::pbrt::*;
use std::sync::Arc;

/// Models scattering properties in volumetric media.
pub trait PhaseFunction {
    /// Returns the value of the phase function for the given pair of directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    fn p(&self, wo: &Vector3f, wi: &Vector3f) -> Float;

    /// Returns the phase function value and sampled incident direction given the
    /// outgoing direction and a sample value in [0, 1)^2.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - Sample value in [0, 1)^2.
    fn sample_p(&self, wo: &Vector3f, u: &Point2f) -> (Float, Vector3f);
}

/// Atomic reference counted `PhaseFunction`.
pub type ArcPhaseFunction = Arc<dyn PhaseFunction + Send + Sync>;
