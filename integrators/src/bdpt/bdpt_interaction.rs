//! BDPT Interactions

/*
#![allow(dead_code)]
use core::geometry::*;
use core::interaction::*;
use core::pbrt::*;
use core::reflection::*;

/// This encapsulates the `core::interaction::Interaction` enums for ease of use here as we don't
/// want to pollute those with `EndpointInteraction`.
///
/// The lifetime specifiers are from `SurfaceInteraction`:
/// * `'scene` - Shared reference to the scene containing primitive.
#[derive(Clone)]
pub(crate) enum BDPTInteraction<'scene> {
    /// Surface interaction.
    Surface { si: SurfaceInteraction<'scene>, bsdf: BSDF },

    /// Medium interaction.
    Medium { mi: MediumInteraction },

    /// Endpoint interaction.
    Endpoint { ei: EndpointInteraction },
}

impl<'scene> Default for BDPTInteraction<'scene> {
    /// Returns the "default value" for `BDPTInteraction<'scene>`.
    fn default() -> Self {
        Self::Endpoint {
            ei: EndpointInteraction::Light {
                hit: Hit::default(),
                light: None,
            },
        }
    }
}

impl<'scene> BDPTInteraction<'scene> {
    /// Spawn's a new ray in the given direction.
    ///
    /// * `d` - The new direction.
    pub(crate) fn spawn_ray(&self, d: &Vector3f) -> Ray {
        match self {
            Self::Surface { si, bsdf: _ } => si.spawn_ray(d),
            Self::Medium { mi } => mi.spawn_ray(d),
            Self::Endpoint { ei } => ei.spawn_ray(d),
        }
    }

    /// Spawn's a new ray towards another point.
    ///
    /// * `p` - The target point.
    pub(crate) fn spawn_ray_to_point(&self, p: &Point3f) -> Ray {
        match self {
            Self::Surface { si, bsdf: _ } => si.spawn_ray_to_point(p),
            Self::Medium { mi } => mi.spawn_ray_to_point(p),
            Self::Endpoint { ei } => ei.spawn_ray_to_point(p),
        }
    }

    /// Spawn's a new ray towards another interaction.
    ///
    /// * `hit` - The interaction.
    pub(crate) fn spawn_ray_to_hit(&self, hit: &Hit) -> Ray {
        match self {
            Self::Surface { si, bsdf: _ } => si.spawn_ray_to_hit(hit),
            Self::Medium { mi } => mi.spawn_ray_to_hit(hit),
            Self::Endpoint { ei } => ei.spawn_ray_to_hit(hit),
        }
    }

    pub(crate) fn p(&self) -> Point3f {
        match self {
            Self::Surface { si, bsdf: _ } => si.hit.p,
            Self::Medium { mi } => mi.hit.p,
            Self::Endpoint { ei } => ei.hit().p,
        }
    }

    pub(crate) fn time(&self) -> Float {
        match self {
            Self::Surface { si, bsdf: _ } => si.hit.time,
            Self::Medium { mi } => mi.hit.time,
            Self::Endpoint { ei } => ei.hit().time,
        }
    }

    pub(crate) fn ng(&self) -> Normal3f {
        match self {
            Self::Surface { si, bsdf: _ } => si.hit.n,
            Self::Medium { mi } => mi.hit.n,
            Self::Endpoint { ei } => ei.hit().n,
        }
    }

    pub(crate) fn ns(&self) -> Normal3f {
        match self {
            Self::Surface { si, bsdf: _ } => si.shading.n,
            Self::Medium { mi } => mi.hit.n,
            Self::Endpoint { ei } => ei.hit().n,
        }
    }
}
*/
