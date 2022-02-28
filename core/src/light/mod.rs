//! Light

#![allow(dead_code)]
use crate::geometry::*;
use crate::interaction::*;
use crate::pbrt::*;
use crate::scene::*;
use crate::spectrum::*;
use std::sync::Arc;

mod light_type;
mod visibility_tester;

/// Return value for `Light::sample_li()`.
#[derive(Clone)]
pub struct Li {
    /// Incident direction.
    pub wi: Vector3f,

    /// PDF.
    pub pdf: Float,

    /// Visibility tester.
    pub visibility: Option<VisibilityTester>,

    /// Radiance arriving at intersection point.
    pub value: Spectrum,
}

impl Li {
    /// Return a new `Li`.
    ///
    /// * `wi`         - Incident direction.
    /// * `pdf`        - PDF.
    /// * `visibility` - Visibility tester.
    /// * `value`      - Radiance arriving at intersection point.
    pub fn new(
        wi: Vector3f,
        pdf: Float,
        visibility: Option<VisibilityTester>,
        value: Spectrum,
    ) -> Self {
        Self {
            wi,
            pdf,
            visibility,
            value,
        }
    }
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

impl Le {
    /// Return a new `Le`.
    ///
    /// * `ray`     - Ray leaving the light source.
    /// * `n_light` - Surface normal at the point on the light source.
    /// * `pdf_pos` - The ray origin's probability density with respect to surface
    ///               area on the light.
    /// * `pdf_dir` - The ray directions's probability density with respect to
    ///               solid angle.
    /// * `value`   - Emitted radiance value.
    pub fn new(
        ray: Ray,
        n_light: Normal3f,
        pdf_pos: Float,
        pdf_dir: Float,
        value: Spectrum,
    ) -> Self {
        Self {
            ray,
            n_light,
            pdf_pos,
            pdf_dir,
            value,
        }
    }
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

impl Pdf {
    /// Return a new `Pdf`.
    ///
    /// * `pdf_pos` - The ray origin's probability density with respect to
    ///               surface area on the light.
    /// * `pdf_dir` - The ray directions's probability density with respect to
    ///               solid angle.
    pub fn new(pdf_pos: Float, pdf_dir: Float) -> Self {
        Self { pdf_pos, pdf_dir }
    }
}

/// Light trait provides common behavior.
pub trait Light {
    /// Initialize the light source before rendering begins.
    ///
    /// * `scene` - The scene.
    fn preprocess(&self, _scene: &Scene) {}

    /// Returns the type of light.
    fn get_type(&self) -> LightType;

    /// Return the radiance arriving at an interaction point.
    ///
    /// * `hit` - The interaction hit point.
    /// * `u`   - Sample value for Monte Carlo integration.
    fn sample_li(&self, hit: &Hit, u: &Point2f) -> Li;

    /// Return the total emitted power.
    fn power(&self) -> Spectrum;

    /// Returns emitted radiance due to that light along a ray that escapes the
    /// scene bounds.
    ///
    /// * `ray` - The ray with differentials.
    fn le(&self, _ray: &Ray) -> Spectrum {
        Spectrum::ZERO
    }

    /// Returns the probability density with respect to solid angle for the light’s
    /// `sample_li()`.
    ///
    /// * `hit` - The interaction hit point.
    /// * `wi`  - The incident direction.
    fn pdf_li(&self, hit: &Hit, wi: &Vector3f) -> Float;

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
    fn is_delta_light(&self) -> bool {
        self.get_type().is_delta_light()
    }

    /// Returns the number of samples to use for the light source.
    fn get_num_samples(&self) -> usize;
}

/// Atomic reference counted `Light`.
pub type ArcLight = Arc<dyn Light + Send + Sync>;

/// AreaLight trait provides common behavior for area lights.
pub trait AreaLight: Light {
    /// Returns the area light's emitted radiance in a given outgoing direction.
    ///
    /// * `it` - Point on a surface to evaluate emitted radiance.
    /// * `w`  - Outgoing direction.
    fn l(&self, hit: &Hit, w: &Vector3f) -> Spectrum;
}

/// Atomic reference counted `AreaLight`.
pub type ArcAreaLight = Arc<dyn AreaLight + Send + Sync>;

// Re-export
pub use light_type::*;
pub use visibility_tester::*;
