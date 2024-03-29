//! Distant Source

use core::geometry::*;
use core::interaction::*;
use core::light::*;
use core::medium::*;
use core::paramset::*;
use core::pbrt::*;
use core::sampling::*;
use core::scene::*;
use core::spectrum::*;
use std::sync::{Arc, RwLock};

/// Implements a directional light source that deposits illumination from the same direction at every point in space.
#[derive(Clone)]
pub struct DistantLight {
    /// Light ID. This is usually the index of the light in the scene's light sources. Usefull for adding lights into
    /// `std::collections::HashMap`.
    pub id: usize,

    /// Light source type.
    pub light_type: LightType,

    /// Participating medium.
    pub medium_interface: MediumInterface,

    /// Transformation from light coordinate system to world coordinate system.
    pub light_to_world: ArcTransform,

    /// The emitted radiance `L`.
    pub emitted_radiance: Spectrum,

    /// Direction of light.
    pub w_light: Vector3f,

    /// Center of the world.
    pub world_center: Arc<RwLock<Point3f>>,

    /// Radius of the spherical world bounds.
    pub world_radius: Arc<RwLock<Float>>,
}

impl DistantLight {
    /// Returns a new `DistantLight`.
    ///
    /// * `id`               - Light ID.
    /// * `light_to_world`   - Transformation from light coordinate system to world coordinate system.
    /// * `emitted_radiance` - The emitted radiance.
    /// * `w_light`          - Direction of light.
    pub fn new(id: usize, light_to_world: ArcTransform, emitted_radiance: Spectrum, w_light: Vector3f) -> Self {
        let light_to_world = Arc::clone(&light_to_world);
        let w_light = light_to_world.transform_vector(&w_light).normalize();

        Self {
            id,
            light_type: LightType::DELTA_DIRECTION_LIGHT,
            light_to_world,
            medium_interface: MediumInterface::vacuum(),
            world_center: Arc::new(RwLock::new(Point3f::ZERO)), // Calculated in preprocess().
            world_radius: Arc::new(RwLock::new(1.0)),           // Calculated in preprocess().
            w_light,
            emitted_radiance,
        }
    }
}

impl Light for DistantLight {
    /// Initialize the light source before rendering begins.
    ///
    /// * `scene` - The scene.
    fn preprocess(&self, scene: &Scene) {
        let (world_center, world_radius) = scene.world_bound.bounding_sphere();
        *self.world_center.write().unwrap() = world_center;
        *self.world_radius.write().unwrap() = world_radius;
    }

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
        let world_radius = *self.world_radius.read().unwrap();
        let p_outside = hit.p + self.w_light * (2.0 * world_radius);

        let p0 = hit.clone();
        let p1 = Hit::new_minimal(p_outside, hit.time, Some(self.medium_interface.clone()));
        let vis = VisibilityTester::new(p0, p1);

        Some(Li::new(self.w_light, 1.0, vis, self.emitted_radiance))
    }

    /// Return the total emitted power.
    fn power(&self) -> Spectrum {
        let world_radius = *self.world_radius.read().unwrap();
        self.emitted_radiance * PI * world_radius * world_radius
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
        let world_center = *self.world_center.read().unwrap();
        let world_radius = *self.world_radius.read().unwrap();

        // Choose point on disk oriented toward infinite light direction.
        let (v1, v2) = coordinate_system(&self.w_light);
        let cd = concentric_sample_disk(u1);
        let p_disk = world_center + world_radius * (cd.x * v1 + cd.y * v2);

        // Set ray origin and direction for infinite light ray.
        let dir = -self.w_light;
        let ray = Ray::new(
            p_disk + world_radius * self.w_light,
            dir,
            INFINITY,
            time,
            self.medium_interface.inside.as_ref().map(Arc::clone),
        );
        Le::new(
            ray,
            Normal3f::from(dir),
            1.0 / (PI * world_radius * world_radius),
            1.0,
            self.emitted_radiance,
        )
    }

    /// Returns the probability density for the light’s `sample_le()`.
    ///
    /// * `ray`     - The ray.
    /// * `n_light` - The normal.
    fn pdf_le(&self, _ray: &Ray, _n_light: &Normal3f) -> Pdf {
        let world_radius = *self.world_radius.read().unwrap();
        Pdf::new(1.0 / (PI * world_radius * world_radius), 0.0)
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

impl From<(&ParamSet, ArcTransform, usize)> for DistantLight {
    /// Create a `DistantLight` from given parameter set, light to world transform and id.
    ///
    /// * `p` - A tuple containing the parameter set, light to world transform and id.
    fn from(p: (&ParamSet, ArcTransform, usize)) -> Self {
        let (params, light_to_world, id) = p;

        let emitted_radiance = params.find_one_spectrum("L", Spectrum::ONE);
        let sc = params.find_one_spectrum("scale", Spectrum::ONE);
        let from = params.find_one_point3f("from", Point3f::new(0.0, 0.0, 0.0));
        let to = params.find_one_point3f("to", Point3f::new(0.0, 0.0, 1.0));
        let dir = from - to;
        Self::new(id, Arc::clone(&light_to_world), emitted_radiance * sc, dir)
    }
}
