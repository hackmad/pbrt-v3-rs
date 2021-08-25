//! Medium Interactions

#![allow(dead_code)]
use crate::geometry::*;
use crate::medium::*;
use crate::pbrt::*;
use std::sync::Arc;

/// MediumInteraction represents an interaction point in a scattering medium.
#[derive(Clone)]
pub struct MediumInteraction {
    /// The common interaction data.
    pub hit: Hit,

    /// The phase function.
    pub phase: ArcPhaseFunction,
}

impl MediumInteraction {
    /// Create a new medium interaction.
    ///
    /// * `p`      - The point of interaction.
    /// * `wo`     - The negative ray direction (outgoing direction used
    /// *            when computing lighting at points).
    /// * `time`   - Time when interaction occurred.
    /// * `medium` - The medium.
    /// * `phase`  - The phase function.
    pub fn new(
        p: Point3f,
        wo: Vector3f,
        time: Float,
        medium: ArcMedium,
        phase: ArcPhaseFunction,
    ) -> Self {
        Self {
            hit: Hit::new(
                p,
                time,
                Vector3f::default(),
                wo,
                Normal3f::default(),
                Some(MediumInterface::from(medium)),
            ),
            phase: Arc::clone(&phase),
        }
    }
}
