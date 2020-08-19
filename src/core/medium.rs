//! Medium

#![allow(dead_code)]
use std::sync::Arc;

/// Medium trait to handle volumetric scattering properties.
pub trait Medium {}

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

/// Create a medium interface between two media.
///
/// * `inside`  - The interior medium.
/// * `outside` - The exterior medium.
pub fn medium_interface(inside: ArcMedium, outside: ArcMedium) -> MediumInterface {
    MediumInterface {
        inside: inside.clone(),
        outside: outside.clone(),
    }
}

/// Create a medium interface between same media.
///
/// * `medium`  - The medium on either side of the interface.
pub fn medium_interface_same(medium: ArcMedium) -> MediumInterface {
    MediumInterface {
        inside: medium.clone(),
        outside: medium.clone(),
    }
}

impl MediumInterface {
    /// Returns `true` if the medium interface marks a transition between
    /// two distinct media.
    pub fn is_medium_transition(&self) -> bool {
        Arc::ptr_eq(&self.inside, &self.outside)
    }
}
