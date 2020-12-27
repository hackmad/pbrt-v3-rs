//! Medium

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::sampler::*;
use crate::core::spectrum::*;
use std::sync::Arc;

/// Medium trait to handle volumetric scattering properties.
pub trait Medium {
    /// Returns the beam transmittance along a given ray.
    ///
    /// * `ray`     - The ray.
    /// * `sampler` - The sampler.
    fn tr(&self, ray: &Ray, sampler: ArcSampler) -> Spectrum;
}

/// Atomic reference counted `Medium`.
pub type ArcMedium = Arc<dyn Medium + Send + Sync>;

/// MediumInterface represents the boundary interface between two media.
#[derive(Clone)]
pub struct MediumInterface {
    /// Represent the interior of a geometric primitive.
    pub inside: ArcMedium,

    /// Represent the exterior of a geometric primitive.
    pub outside: ArcMedium,
}

impl MediumInterface {
    /// Create a medium interface between two media.
    ///
    /// * `inside`  - The interior medium.
    /// * `outside` - The exterior medium.
    pub fn new(inside: ArcMedium, outside: ArcMedium) -> Self {
        Self {
            inside: inside.clone(),
            outside: outside.clone(),
        }
    }

    /// Returns `true` if the medium interface marks a transition between
    /// two distinct media.
    pub fn is_medium_transition(&self) -> bool {
        Arc::ptr_eq(&self.inside, &self.outside)
    }
}

impl From<ArcMedium> for MediumInterface {
    /// Create a medium interface between same media.
    ///
    /// * `medium` - The medium on either side of the interface.
    fn from(medium: ArcMedium) -> Self {
        Self {
            inside: medium.clone(),
            outside: medium.clone(),
        }
    }
}
