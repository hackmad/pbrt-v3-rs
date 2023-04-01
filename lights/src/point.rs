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

/// Implements an isotropic point light source that emits the same amount of light in all directions.
#[derive(Clone)]
pub struct PointLight {
    /// Light ID. This is usually the index of the light in the scene's light sources. Usefull for adding lights into
    /// `std::collections::HashMap`.
    pub id: usize,

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
    /// * `id`               - Light ID.
    /// * `light_to_world`   - Transformation from light coordinate system to world coordinate system.
    /// * `medium_interface` - Participating medium.
    /// * `intensity`        - Intensity.
    pub fn new(
        id: usize,
        light_to_world: ArcTransform,
        medium_interface: MediumInterface,
        intensity: Spectrum,
    ) -> Self {
        let light_to_world = Arc::clone(&light_to_world);
        let world_to_light = Arc::new(light_to_world.inverse());
        let p_light = light_to_world.transform_point(&Point3f::ZERO);

        Self {
            id,
            light_type: LightType::DELTA_POSITION_LIGHT,
            medium_interface: medium_interface.clone(),
            light_to_world,
            world_to_light,
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

    /// Returns the light unique id. Usually the index in the scene's light sources.
    fn get_id(&self) -> usize {
        self.id
    }

    /// Return the radiance arriving at an interaction point.
    ///
    /// * `hit` - The interaction hit point.
    /// * `u`   - Sample value for Monte Carlo integration.
    fn sample_li(&self, hit: &Hit, _u: &Point2f) -> Option<Li> {
        let wi = (self.p_light - hit.p).normalize();
        let pdf = 1.0;

        let p0 = hit.clone();
        let p1 = Hit::new_minimal(self.p_light, hit.time, Some(self.medium_interface.clone()));
        let vis = VisibilityTester::new(p0, p1);

        let value = self.intensity / self.p_light.distance_squared(hit.p);
        Some(Li::new(wi, pdf, vis, value))
    }

    /// Return the total emitted power.
    fn power(&self) -> Spectrum {
        FOUR_PI * self.intensity
    }

    /// Returns the probability density with respect to solid angle for the light’s `sample_li()`.
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
            self.medium_interface.inside.as_ref().map(Arc::clone),
        );
        Le::new(ray, Normal3f::from(dir), 1.0, uniform_sphere_pdf(), self.intensity)
    }

    /// Returns the probability density for the light’s `sample_le()`.
    ///
    /// * `ray`     - The ray.
    /// * `n_light` - The normal.
    fn pdf_le(&self, _ray: &Ray, _n_light: &Normal3f) -> Pdf {
        Pdf::new(0.0, uniform_sphere_pdf())
    }

    /// Returns the number of samples to use for the light source.
    fn get_num_samples(&self) -> usize {
        1
    }

    /// Returns the area light's emitted radiance in a given outgoing direction.
    ///
    /// * `it` - Point on a surface to evaluate emitted radiance.
    /// * `w`  - Outgoing direction.
    fn l(&self, _hit: &Hit, _w: &Vector3f) -> Spectrum {
        panic!("Invalid call to Light::l() for non area lights.")
    }
}

impl From<(&ParamSet, ArcTransform, Option<ArcMedium>, usize)> for PointLight {
    /// Create a `PointLight` from given parameter set, light to world transform, medium and id.
    ///
    /// * `p` - A tuple containing the parameter set, light to world transform, medium and id.
    fn from(p: (&ParamSet, ArcTransform, Option<ArcMedium>, usize)) -> Self {
        let (params, light_to_world, medium, id) = p;

        let intensity = params.find_one_spectrum("I", Spectrum::ONE);
        let sc = params.find_one_spectrum("scale", Spectrum::ONE);
        let p = params.find_one_point3f("from", Point3f::ZERO);
        let l2w = Transform::translate(&Vector3f::new(p.x, p.y, p.z)) * light_to_world.as_ref();
        Self::new(id, Arc::new(l2w), MediumInterface::from(medium), intensity * sc)
    }
}
