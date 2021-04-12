//! Light

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::scene::*;
use crate::core::spectrum::*;
use std::sync::Arc;

mod light_type;
mod visibility_tester;

/// Return value for `Light::sample_li()`.
#[derive(Clone)]
pub struct Li {
    /// Incoming direction.
    pub wi: Vector3f,

    /// PDF.
    pub pdf: Float,

    /// Visibility tester.
    pub visibility: VisibilityTester,

    /// Radiance arriving at intersection point.
    pub value: Spectrum,
}

/// Return value for `Light::sample_le()`.
#[derive(Clone)]
pub struct Le {
    /// Ray leaving the light source.
    pub ray: Ray,

    /// Surface normal at the point on the light source.
    pub n_light: Normal3f,

    /// The ray origin's probability density with respect to surface area on the
    /// light.
    pub pdf_pos: Float,

    /// The ray directions's probability density with respect to solid angle.
    pub pdf_dir: Float,

    /// Emitted radiance value.
    pub value: Spectrum,
}

/// Return value for `Light::pdf_le()`.
#[derive(Copy, Clone)]
pub struct Pdf {
    /// The ray origin's probability density with respect to surface area on the
    /// light.
    pub pdf_pos: Float,

    /// The ray directions's probability density with respect to solid angle.
    pub pdf_dir: Float,
}

/// Light trait provides common behavior.
pub trait Light {
    /// Initialize the light source before rendering begins.
    ///
    /// * `scene` - The scene.
    fn preprocess(&mut self, _scene: &Scene) {}

    /// Returns the type of light.
    fn get_type(&self) -> LightType;

    /// Return the radiance arriving at an interaction point.
    ///
    /// * `reference` - The interaction point.
    /// * `u`         - Sample value for Monte Carlo integration.
    fn sample_li(&self, reference: &dyn Interaction, u: &Point2f) -> Li;

    /// Return the total emitted power.
    fn power(&self) -> Spectrum;

    /// Returns emitted radiance due to that light along a ray that escapes the
    /// scene bounds.
    ///
    /// * `r` - The ray differentials.
    fn le(&self, r: &RayDifferential) -> Spectrum;

    /// Returns the probability density with respect to solid angle for the light’s
    /// `sample_li()`.
    ///
    /// * `reference` - The interaction point.
    /// * `wi`        - The incoming direction.
    fn pdf_li(&self, reference: ArcInteraction, wi: &Vector3f) -> Float;

    /// Returns a sampled light-carrying ray leaving the light source.
    ///
    /// * `u1`   - Sample values for Monte Carlo.
    /// * `u2`   - Sample values for Monte Carlo.
    /// * `time` - Time to use for the ray.
    fn sample_le(&self, u1: &Point2f, u2: &Point2f, time: Float) -> Le;

    /// Returns the probability density for the light’s `sample_le()`.
    ///
    /// * `ray`     - The ray.
    /// * `n_light` - The normal.
    fn pdf_le(&self, ray: &Ray, n_light: &Normal3f) -> Pdf;

    /// Returns whether light source is a delta light.
    fn is_delta_light(&self) -> bool;
}

/// Atomic reference counted `Light`.
pub type ArcLight = Arc<dyn Light + Send + Sync>;

/// AreaLight trait provides common behavior for area lights.
pub trait AreaLight: Light {
    /// Return struct as `Light`.
    fn as_light(&self) -> &'static dyn Light;

    /// Returns the area light's emitted radiance in a given outgoing direction.
    ///
    /// * `it` - Point on a surface to evaluate emitted radiance.
    /// * `w`  - Outgoing direction.
    fn l(&self, it: &dyn Interaction, w: &Vector3f) -> Spectrum;
}

/// Atomic reference counted `AreaLight`.
pub type ArcAreaLight = Arc<dyn AreaLight + Send + Sync>;

// Re-export
pub use light_type::*;
pub use visibility_tester::*;
