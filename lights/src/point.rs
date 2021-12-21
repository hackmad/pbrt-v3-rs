//! Point Light Source

use core::geometry::*;
use core::interaction::*;
use core::light::*;
use core::medium::*;
use core::paramset::*;
use core::pbrt::*;
use core::sampling::*;
use core::spectrum::*;
use std::sync::Arc;

/// Implements an isotropic point light source that emits the same amount of
/// light in all directions.
#[derive(Clone)]
pub struct PointLight {
    /// Light source type.
    pub light_type: LightType,

    /// Participating medium.
    pub medium_interface: MediumInterface,

    /// Transformation from light coordinate system to world coordinate system.
    pub light_to_world: ArcTransform,

    /// Transformation from world coordinate system to light coordinate system.
    pub world_to_light: ArcTransform,

    /// Position.
    pub p_light: Point3f,

    /// Intensity.
    pub intensity: Spectrum,
}

impl PointLight {
    /// Returns a new `PointLight`.
    ///
    /// * `light_to_world`   - Transformation from light coordinate system to
    ///                        world coordinate system.
    /// * `medium_interface` - Participating medium.
    /// * `intensity`        - Intensity.
    pub fn new(
        light_to_world: ArcTransform,
        medium_interface: MediumInterface,
        intensity: Spectrum,
    ) -> Self {
        let world_to_light = Arc::clone(&light_to_world).inverse();
        let p_light = Arc::clone(&light_to_world).transform_point(&Point3f::default());
        Self {
            light_type: LightType::DELTA_POSITION_LIGHT,
            medium_interface: medium_interface.clone(),
            light_to_world: Arc::clone(&light_to_world),
            world_to_light: Arc::new(world_to_light),
            p_light,
            intensity,
        }
    }
}

impl Light for PointLight {
    /// Returns the type of light.
    fn get_type(&self) -> LightType {
        self.light_type
    }

    /// Return the radiance arriving at an interaction point.
    ///
    /// * `hit` - The interaction hit point.
    /// * `u`   - Sample value for Monte Carlo integration.
    fn sample_li(&self, hit: &Hit, _u: &Point2f) -> Li {
        let wi = (self.p_light - hit.p).normalize();
        let pdf = 1.0;

        let p0 = hit.clone();
        let p1 = Hit::new_minimal(self.p_light, hit.time, hit.medium_interface.clone());
        let vis = VisibilityTester::new(p0, p1);

        let value = self.intensity / self.p_light.distance_squared(hit.p);
        Li::new(wi, pdf, Some(vis), value)
    }

    /// Return the total emitted power.
    fn power(&self) -> Spectrum {
        FOUR_PI * self.intensity
    }

    /// Returns the probability density with respect to solid angle for the light’s
    /// `sample_li()`.
    ///
    /// * `hit` - The interaction hit point.
    /// * `wi`  - The incident direction.
    fn pdf_li(&self, _hit: &Hit, _wi: &Vector3f) -> Float {
        0.0
    }

    /// Returns a sampled light-carrying ray leaving the light source.
    ///
    /// * `u1`   - Sample values for Monte Carlo.
    /// * `u2`   - Sample values for Monte Carlo.
    /// * `time` - Time to use for the ray.
    fn sample_le(&self, u1: &Point2f, _u2: &Point2f, time: Float) -> Le {
        let dir = uniform_sample_sphere(&u1);
        let ray = Ray::new(
            self.p_light,
            dir,
            INFINITY,
            time,
            self.medium_interface.inside.clone(),
        );
        Le::new(
            ray,
            Normal3f::from(dir),
            1.0,
            uniform_sphere_pdf(),
            self.intensity,
        )
    }

    /// Returns the probability density for the light’s `sample_le()`.
    ///
    /// * `ray`     - The ray.
    /// * `n_light` - The normal.
    fn pdf_le(&self, _ray: &Ray, _n_light: &Normal3f) -> Pdf {
        Pdf::new(0.0, uniform_sphere_pdf())
    }
}

impl From<(&ParamSet, ArcTransform, Option<ArcMedium>)> for PointLight {
    /// Create a `PointLight` from given parameter set, light to world transform
    /// and medium.
    ///
    /// * `p` - A tuple containing the parameter set, light to world transform
    ///         and medium.
    fn from(p: (&ParamSet, ArcTransform, Option<ArcMedium>)) -> Self {
        let (params, light_to_world, medium) = p;

        let intensity = params.find_one_spectrum("I", Spectrum::new(1.0));
        let sc = params.find_one_spectrum("scale", Spectrum::new(1.0));
        let p = params.find_one_point3f("from", Point3f::default());
        let l2w = Transform::translate(&Vector3f::new(p.x, p.y, p.z)) * light_to_world.as_ref();
        Self::new(Arc::new(l2w), MediumInterface::from(medium), intensity * sc)
    }
}
