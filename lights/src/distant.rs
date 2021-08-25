//! Distant Source

use core::geometry::*;
use core::light::*;
use core::medium::*;
use core::paramset::*;
use core::pbrt::*;
use core::sampling::*;
use core::scene::*;
use core::spectrum::*;
use std::sync::Arc;

/// Implements a directional light source that deposits illumination from the
/// same direction at every point in space.
#[derive(Clone)]
pub struct DistantLight {
    /// Light source type.
    pub light_type: LightType,

    /// Participating medium.
    pub medium_interface: MediumInterface,

    /// Transformation from light coordinate system to world coordinate system.
    pub light_to_world: ArcTransform,

    /// Transformation from world coordinate system to light coordinate system.
    pub world_to_light: ArcTransform,

    /// The emitted radiance `L`.
    pub emitted_radiance: Spectrum,

    /// Direction of light.
    pub w_light: Vector3f,

    /// Center of the world.
    pub world_center: Point3f,

    /// Radius of the spherical world bounds.
    pub world_radius: Float,
}

impl DistantLight {
    /// Returns a new `DistantLight`.
    ///
    /// * `light_to_world`   - Transformation from light coordinate system to
    ///                        world coordinate system.
    /// * `emitted_radiance` - The emitted radiance.
    /// * `w_light`          - Direction of light.
    pub fn new(
        light_to_world: ArcTransform,
        emitted_radiance: Spectrum,
        w_light: Vector3f,
    ) -> Self {
        let world_to_light = Arc::clone(&light_to_world).inverse();

        Self {
            light_type: LightType::from(DELTA_DIRECTION_LIGHT),
            light_to_world: Arc::clone(&light_to_world),
            world_to_light: Arc::new(world_to_light),
            medium_interface: MediumInterface::vacuum(),
            world_center: Point3f::default(), // Calculated in preprocess().
            world_radius: 1.0,                // Calculated in preprocess().
            w_light,
            emitted_radiance,
        }
    }
}

impl Light for DistantLight {
    /// Initialize the light source before rendering begins.
    ///
    /// * `scene` - The scene.
    fn preprocess(&mut self, scene: &Scene) {
        let (world_center, world_radius) = scene.world_bound.bounding_sphere();
        self.world_center = world_center;
        self.world_radius = world_radius;
    }

    /// Returns the type of light.
    fn get_type(&self) -> LightType {
        self.light_type
    }

    /// Return the radiance arriving at an interaction point.
    ///
    /// * `hit` - The interaction hit point.
    /// * `u`   - Sample value for Monte Carlo integration.
    fn sample_li(&self, hit: &Hit, _u: &Point2f) -> Li {
        let wi = self.w_light;
        let pdf = 1.0;
        let p_outside = hit.p + self.w_light * (2.0 * self.world_radius);
        let visibility = Some(VisibilityTester::new(hit.clone(), p_outside));
        let value = self.emitted_radiance;
        Li::new(wi, pdf, visibility, value)
    }

    /// Return the total emitted power.
    fn power(&self) -> Spectrum {
        self.emitted_radiance * PI * self.world_radius * self.world_radius
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
        // Choose point on disk oriented toward infinite light direction.
        let (v1, v2) = coordinate_system(&self.w_light);
        let cd = concentric_sample_disk(u1);
        let p_disk = self.world_center + self.world_radius * (cd.x * v1 + cd.y * v2);

        // Set ray origin and direction for infinite light ray.
        let dir = -self.w_light;
        let ray = Ray::new(
            p_disk + self.world_radius * self.w_light,
            dir,
            INFINITY,
            time,
            self.medium_interface.inside.clone(),
        );
        Le::new(
            ray,
            Normal3f::from(dir),
            1.0 / (PI * self.world_radius * self.world_radius),
            1.0,
            self.emitted_radiance,
        )
    }

    /// Returns the probability density for the light’s `sample_le()`.
    ///
    /// * `ray`     - The ray.
    /// * `n_light` - The normal.
    fn pdf_le(&self, _ray: &Ray, _n_light: &Normal3f) -> Pdf {
        Pdf::new(1.0 / (PI * self.world_radius * self.world_radius), 0.0)
    }
}

impl From<(&ParamSet, ArcTransform)> for DistantLight {
    /// Create a `DistantLight` from given parameter set and light to world transform.
    ///
    /// * `p` - A tuple containing the parameter set and light to world transform.
    fn from(p: (&ParamSet, ArcTransform)) -> Self {
        let (params, light_to_world) = p;

        let emitted_radiance = params.find_one_spectrum("L", Spectrum::new(1.0));
        let sc = params.find_one_spectrum("scale", Spectrum::new(1.0));
        let from = params.find_one_point3f("from", Point3f::new(0.0, 0.0, 0.0));
        let to = params.find_one_point3f("to", Point3f::new(0.0, 0.0, 0.1));
        let dir = from - to;
        Self::new(Arc::clone(&light_to_world), emitted_radiance * sc, dir)
    }
}
