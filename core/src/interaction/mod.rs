//! Interactions

#![allow(dead_code)]
use crate::geometry::*;
use crate::medium::*;
use crate::pbrt::*;
use std::fmt;
use std::sync::Arc;

mod medium_interaction;
mod surface_interaction;

pub use medium_interaction::*;
use num_traits::Zero;
pub use surface_interaction::*;

/// Interaction enumeration.
///
/// The lifetime specifiers are from `SurfaceInteraction`:
/// * `'scene` - Shared reference to the scene containing primitive.
/// * `'arena` - Shared mutable reference to values allocated by a memory arena.
///              Because we use bumpalo::Bump, this is how it returns allocated
///              values.
pub enum Interaction<'scene, 'arena> {
    /// Represents geometry of a particular point on a surface.
    Surface { si: SurfaceInteraction<'scene> },

    /// Represents an interaction point in a scattering medium.
    Medium { mi: MediumInteraction<'arena> },
}

impl<'scene, 'arena> Interaction<'scene, 'arena> {
    /// Returns the interaction hit point.
    pub fn get_hit(&self) -> &Hit {
        match self {
            Self::Surface { si } => &si.hit,
            Self::Medium { mi } => &mi.hit,
        }
    }

    /// Spawn's a new ray in the given direction.
    ///
    /// * `d` - The new direction.
    pub fn spawn_ray(&self, d: &Vector3f) -> Ray {
        match self {
            Self::Surface { si } => si.spawn_ray(d),
            Self::Medium { mi } => mi.spawn_ray(d),
        }
    }

    /// Spawn's a new ray towards another point.
    ///
    /// * `p` - The target point.
    pub fn spawn_ray_to_point(&self, p: &Point3f) -> Ray {
        match self {
            Self::Surface { si } => si.spawn_ray_to_point(p),
            Self::Medium { mi } => mi.spawn_ray_to_point(p),
        }
    }

    /// Spawn's a new ray towards another interaction.
    ///
    /// * `hit` - The interaction.
    pub fn spawn_ray_to_hit(&self, hit: &Hit) -> Ray {
        match self {
            Self::Surface { si } => si.spawn_ray_to_hit(hit),
            Self::Medium { mi } => mi.spawn_ray_to_hit(hit),
        }
    }
}

/// Hit provides common data shared by implementations of `Interaction` trait.
#[derive(Clone)]
pub struct Hit {
    /// Point of interaction.
    pub p: Point3f,

    /// Time when interaction occurred.
    pub time: Float,

    /// Floating point error for ray intersection points.
    pub p_error: Vector3f,

    /// The negative ray direction (outgoing direction used when computing
    /// lighting at points).
    pub wo: Vector3f,

    /// Surface normal at the point `p`.
    pub n: Normal3f,

    /// The medium interface used for scattering media.
    pub medium_interface: Option<MediumInterface>,
}

/// Atomic reference counted `Hit`.
type ArcHit = Arc<Hit>;

impl Hit {
    /// Create a new hit.
    ///
    /// NOTE: If you need to contruct a new `Hit` without `wo`, `n` and `p_error`
    /// use `Hit::new_minimal()`. This function calls `wo.normalize()` and will
    /// generate weird values for zero vectors.
    ///
    /// `p`                - Point of interaction.
    /// `time`             - Time when interaction occurred.
    /// `p_error`          - Floating point error for ray intersection points.
    /// `wo`               - The negative ray direction (outgoing direction used
    ///                      when computing lighting at points).
    /// `n`                - Surface normal at the point `p`.
    /// `medium_interface` - The medium interface used for scattering media.
    pub fn new(
        p: Point3f,
        time: Float,
        p_error: Vector3f,
        wo: Vector3f,
        n: Normal3f,
        medium_interface: Option<MediumInterface>,
    ) -> Self {
        let l2 = wo.length_squared();
        let wo = if l2.is_zero() { wo } else { wo / l2.sqrt() };

        Self {
            p,
            time,
            p_error,
            wo,
            n,
            medium_interface,
        }
    }

    /// Create a new hit from minimal fields.
    ///
    /// `p`                - Point of interaction.
    /// `time`             - Time when interaction occurred.
    /// `medium_interface` - The medium interface used for scattering media.
    pub fn new_minimal(p: Point3f, time: Float, medium_interface: Option<MediumInterface>) -> Self {
        Self {
            p,
            time,
            p_error: Vector3f::ZERO,
            wo: Vector3f::ZERO,
            n: Normal3f::ZERO,
            medium_interface,
        }
    }

    /// Returns `true` if this is a surface interaction.
    pub fn is_surface_interaction(&self) -> bool {
        self.n != Normal3f::zero()
    }

    /// Returns `true` if this is a medium interaction.
    pub fn is_medium_interaction(&self) -> bool {
        !self.is_surface_interaction()
    }

    /// Spawn's a new ray in the given direction.
    ///
    /// * `d` - The new direction.
    pub fn spawn_ray(&self, d: &Vector3f) -> Ray {
        let origin = Ray::offset_origin(&self.p, &self.p_error, &self.n, d);
        Ray::new(
            origin,
            *d,
            INFINITY,
            self.time,
            self.get_medium_in_direction(d),
        )
    }

    /// Spawn's a new ray towards another point.
    ///
    /// * `p` - The target point.
    pub fn spawn_ray_to_point(&self, p: &Point3f) -> Ray {
        let d = p - self.p;
        let origin = Ray::offset_origin(&self.p, &self.p_error, &self.n, &d);
        Ray::new(
            origin,
            d,
            1.0 - SHADOW_EPSILON,
            self.time,
            self.get_medium_in_direction(&d),
        )
    }

    /// Spawn's a new ray towards another interaction.
    ///
    /// * `hit` - The interaction.
    pub fn spawn_ray_to_hit(&self, hit: &Hit) -> Ray {
        let origin = Ray::offset_origin(&self.p, &self.p_error, &self.n, &(hit.p - self.p));
        let target = Ray::offset_origin(&hit.p, &hit.p_error, &hit.n, &(origin - hit.p));
        let d = target - origin;
        Ray::new(
            origin,
            d,
            1.0 - SHADOW_EPSILON,
            self.time,
            self.get_medium_in_direction(&d),
        )
    }

    /// Returns the medium towards a direction.
    ///
    /// * `w` - The direction.
    pub fn get_medium_in_direction(&self, w: &Vector3f) -> Option<ArcMedium> {
        if let Some(mi) = self.medium_interface.as_ref() {
            if w.dot(&self.n) > 0.0 {
                mi.outside.as_ref().map(Arc::clone)
            } else {
                mi.inside.as_ref().map(Arc::clone)
            }
        } else {
            None
        }
    }

    /// Returns the medium when interior and exterior are the same.
    pub fn get_medium(&self) -> Option<ArcMedium> {
        if let Some(mi) = self.medium_interface.as_ref() {
            if mi.is_medium_transition() {
                mi.inside.as_ref().map(Arc::clone)
            } else {
                None
            }
        } else {
            None
        }
    }
}

impl fmt::Display for Hit {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Hit {{ p: {}, time: {}, p_error: {}, wo: {}, n: {}, medium_interface: ... }}",
            self.p, self.time, self.p_error, self.wo, self.n,
        )
    }
}
