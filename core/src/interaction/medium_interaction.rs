//! Medium Interactions

#![allow(dead_code)]
use super::Hit;
use crate::geometry::*;
use crate::medium::*;
use crate::pbrt::*;

/// MediumInteraction represents an interaction point in a scattering medium.
pub struct MediumInteraction<'arena> {
    /// The common interaction data.
    pub hit: Hit,

    /// The phase function.
    pub phase: &'arena mut PhaseFunction<'arena>,
}

impl<'arena> MediumInteraction<'arena> {
    /// Create a new medium interaction.
    ///
    /// * `p`      - The point of interaction.
    /// * `wo`     - The negative ray direction (outgoing direction used
    ///              hen computing lighting at points).
    /// * `time`   - Time when interaction occurred.
    /// * `medium` - The medium.
    /// * `phase`  - The phase function.
    pub fn new(
        p: Point3f,
        wo: Vector3f,
        time: Float,
        medium: Option<ArcMedium>,
        phase: &'arena mut PhaseFunction<'arena>,
    ) -> Self {
        Self {
            hit: Hit::new(
                p,
                time,
                Vector3f::ZERO,
                wo,
                Normal3f::ZERO,
                medium.map(MediumInterface::from),
            ),
            phase,
        }
    }
}
