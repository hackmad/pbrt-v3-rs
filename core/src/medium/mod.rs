//! Medium

#![allow(dead_code)]
use crate::geometry::*;
use crate::interaction::MediumInteraction;
use crate::sampler::*;
use crate::spectrum::*;
use bumpalo::Bump;
use std::sync::Arc;

mod henyey_greenstein;
mod measured_ss;
mod phase_function;

// Re-exports
pub use henyey_greenstein::*;
pub use measured_ss::*;
pub use phase_function::*;

/// Medium trait to handle volumetric scattering properties.
pub trait Medium {
    /// Returns the beam transmittance along a given ray.
    ///
    /// * `ray`     - The ray.
    /// * `sampler` - The sampler.
    fn tr(&self, ray: &Ray, sampler: &mut ArcSampler) -> Spectrum;

    /// Samples a medium scattering interaction along a world-space ray.
    ///
    /// The ray will generally have been intersected against the scene geometry;
    /// thus, implementations of this method shouldn’t ever sample a medium
    /// interaction at a point on the ray beyond its `t_max` value.
    ///
    /// NOTE: Calling code will need to assign this medium as we cannot pass
    /// back and `ArcMedium` out of here for `Self`.
    ///
    /// * `arena`   - The arena for memory allocations.
    /// * `ray`     - The ray.
    /// * `sampler` - The sampler.
    fn sample<'arena>(
        &self,
        arena: &'arena Bump,
        ray: &Ray,
        sampler: &mut ArcSampler,
    ) -> (Spectrum, Option<MediumInteraction<'arena>>);
}

/// Atomic reference counted `Medium`.
pub type ArcMedium = Arc<dyn Medium + Send + Sync>;

/// MediumInterface represents the boundary interface between two media.
pub struct MediumInterface {
    /// Represent the interior of a geometric primitive.
    pub inside: Option<ArcMedium>,

    /// Represent the exterior of a geometric primitive.
    pub outside: Option<ArcMedium>,
}

impl Clone for MediumInterface {
    /// Returns a copy of the value.
    fn clone(&self) -> Self {
        Self {
            inside: self.inside.as_ref().map(Arc::clone),
            outside: self.outside.as_ref().map(Arc::clone),
        }
    }
}

impl MediumInterface {
    /// Create a medium interface between two media.
    ///
    /// * `inside`  - The interior medium.
    /// * `outside` - The exterior medium.
    pub fn new(inside: Option<ArcMedium>, outside: Option<ArcMedium>) -> Self {
        Self {
            inside: inside.as_ref().map(Arc::clone),
            outside: outside.as_ref().map(Arc::clone),
        }
    }

    /// Create a medium interface that represents a vacuum.
    pub fn vacuum() -> Self {
        Self {
            inside: None,
            outside: None,
        }
    }

    /// Returns `true` if the medium interface marks a transition between
    /// two distinct media.
    pub fn is_medium_transition(&self) -> bool {
        match (self.inside.clone(), self.outside.clone()) {
            (Some(inside), Some(outside)) => Arc::ptr_eq(&inside, &outside),
            (Some(_), None) => true,
            (None, Some(_)) => true,
            (None, None) => false,
        }
    }
}

impl From<ArcMedium> for MediumInterface {
    /// Create a medium interface between same media.
    ///
    /// * `medium` - The medium on either side of the interface.
    fn from(medium: ArcMedium) -> Self {
        Self {
            inside: Some(Arc::clone(&medium)),
            outside: Some(Arc::clone(&medium)),
        }
    }
}

impl From<Option<ArcMedium>> for MediumInterface {
    /// Create a medium interface between same media.
    ///
    /// * `medium` - The medium on either side of the interface.
    fn from(medium: Option<ArcMedium>) -> Self {
        Self {
            inside: medium.as_ref().map(Arc::clone),
            outside: medium.as_ref().map(Arc::clone),
        }
    }
}
