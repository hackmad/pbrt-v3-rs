//! Medium Interactions

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::medium::*;
use crate::core::pbrt::*;

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
            phase: phase.clone(),
        }
    }
}

impl Interaction for MediumInteraction {
    /// Returns the interaction hit data.
    fn get_hit(&self) -> &Hit {
        &self.hit
    }

    /// Returns the surface interaction.
    ///
    /// NOTE: This is a hack because I don't want to write a function
    /// for retrieving every field in `MediumInteraction`. If there is a clean
    /// way of retrieving a struct that implements an interface I will get rid
    /// of this.
    fn get_surface_interaction(&self) -> Option<&SurfaceInteraction> {
        None
    }

    /// Returns the medium interaction or None.
    ///
    /// NOTE: This is a hack because I don't want to write a function
    /// for retrieving every field in `MediumInteraction`. If there is a clean
    /// way of retrieving a struct that implements an interface I will get rid
    /// of this.
    fn get_medium_interaction(&self) -> Option<&MediumInteraction> {
        Some(self)
    }
}
