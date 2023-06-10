//! Diffuse Area Light Source

use core::app::options;
use core::geometry::*;
use core::interaction::*;
use core::light::*;
use core::medium::*;
use core::paramset::*;
use core::pbrt::*;
use core::rng::ONE_MINUS_EPSILON;
use core::sampling::*;
use core::spectrum::*;
use std::sync::Arc;

/// Implements a basic area light source with uniform spatial and directional radiance distribution.
#[derive(Clone)]
pub struct DiffuseAreaLight {
    /// Light ID. This is usually the index of the light in the scene's light sources. Usefull for adding lights into
    /// `std::collections::HashMap`.
    pub id: usize,

    /// Light source type.
    pub light_type: LightType,

    /// Used to trace multiple shadow rays to the light to compute soft shadows.
    pub n_samples: usize,

    /// Participating medium.
    pub medium_interface: MediumInterface,

    /// Transformation from light coordinate system to world coordinate system.
    pub light_to_world: ArcTransform,

    /// Transformation from world coordinate system to light coordinate system.
    pub world_to_light: ArcTransform,

    /// Emitted radiance.
    pub l_emit: Spectrum,

    /// Shape describing surface of the light source.
    pub shape: ArcShape,

    /// Surface area of the shape.
    pub area: Float,

    /// Indicates whether light source 2-sided.
    pub two_sided: bool,
}

impl DiffuseAreaLight {
    /// Returns a new `DiffuseAreaLight`.
    ///
    /// * `id`               - Light ID.
    /// * `light_to_world`   - Transformation from light coordinate system to world coordinate system.
    /// * `medium_interface` - Participating medium.
    /// * `l_emit`           - Emitted radiance.
    /// * `n_samples`        - Used to trace multiple shadow rays to the light to compute soft shadows. Default to 1.
    /// * `shape`            - Shape describing surface of the light source.
    /// * `two_sided`        - Indicates whether light source 2-sided.
    pub fn new(
        id: usize,
        light_to_world: ArcTransform,
        medium_interface: MediumInterface,
        l_emit: Spectrum,
        n_samples: usize,
        shape: ArcShape,
        two_sided: bool,
    ) -> Self {
        let light_to_world = Arc::clone(&light_to_world);
        let world_to_light = Arc::new(light_to_world.inverse());
        let area = shape.area();

        // Warn if light has transformation with non-uniform scale, though not for Triangles, since this doesn't matter
        // for them.
        if world_to_light.has_scale() && shape.get_type() != "triangle" {
            warn!(
                "Scaling detected in world to light transformation! 
                The system has numerous assumptions, implicit and explicit, 
                that this transform will have no scale factors in it. 
                Proceed at your own risk; your image may have errors."
            );
        }

        Self {
            id,
            light_type: LightType::AREA_LIGHT,
            medium_interface,
            light_to_world,
            world_to_light,
            l_emit,
            n_samples,
            shape,
            two_sided,
            area,
        }
    }
}

impl Light for DiffuseAreaLight {
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
    fn sample_li(&self, hit: &Hit, u: &Point2f) -> Option<Li> {
        let (mut p_shape, pdf) = self.shape.sample_solid_angle(hit, u);
        p_shape.medium_interface = Some(self.medium_interface.clone());

        let mut wi = p_shape.p - hit.p;
        let wi_len_sq = wi.length_squared();
        if pdf == 0.0 || wi_len_sq == 0.0 {
            return None;
        }
        wi /= wi_len_sq.sqrt(); // Normalize

        let value = self.l(&p_shape, &-wi);
        let vis = VisibilityTester::new(hit.clone(), p_shape);

        Some(Li::new(wi, pdf, vis, value))
    }

    /// Return the total emitted power.
    fn power(&self) -> Spectrum {
        let s = if self.two_sided { 2.0 } else { 1.0 };
        s * self.l_emit * self.area * PI
    }

    /// Returns the probability density with respect to solid angle for the light’s `sample_li()`.
    ///
    /// * `hit` - The interaction hit point.
    /// * `wi`  - The incident direction.
    fn pdf_li(&self, hit: &Hit, wi: &Vector3f) -> Float {
        self.shape.pdf_solid_angle(hit, wi)
    }

    /// Returns a sampled light-carrying ray leaving the light source.
    ///
    /// * `u1`   - Sample values for Monte Carlo.
    /// * `u2`   - Sample values for Monte Carlo.
    /// * `time` - Time to use for the ray.
    fn sample_le(&self, u1: &Point2f, u2: &Point2f, time: Float) -> Le {
        let (mut p_shape_hit, pdf_pos) = self.shape.sample(u1);
        p_shape_hit.medium_interface = Some(self.medium_interface.clone());
        let n_light = p_shape_hit.n;
        let pdf_dir: Float;
        let mut w: Vector3f;

        // Sample a cosine-weighted outgoing direction `w` for area light.
        if self.two_sided {
            let mut u = *u2;
            // Choose a side to sample and then remap u[0] to [0,1] before applying cosine-weighted hemisphere sampling
            // for the chosen side.
            if u[0] < 0.5 {
                u[0] = min(u[0] * 2.0, ONE_MINUS_EPSILON);
                w = cosine_sample_hemisphere(&u);
            } else {
                u[0] = min((u[0] - 0.5) * 2.0, ONE_MINUS_EPSILON);
                w = cosine_sample_hemisphere(&u);
                w.z *= -1.0;
            }
            pdf_dir = 0.5 * cosine_hemisphere_pdf(abs(w.z));
        } else {
            w = cosine_sample_hemisphere(&u2);
            pdf_dir = cosine_hemisphere_pdf(w.z);
        }

        let n = Vector3f::from(p_shape_hit.n);
        let (v1, v2) = coordinate_system(&n);
        w = w.x * v1 + w.y * v2 + w.z * n;

        let mut ray = p_shape_hit.spawn_ray(&w);
        ray.time = time;

        let value = self.l(&p_shape_hit, &w);
        Le::new(ray, n_light, pdf_pos, pdf_dir, value)
    }

    /// Returns the probability density for the light’s `sample_le()`.
    ///
    /// * `ray`     - The ray.
    /// * `n_light` - The normal.
    fn pdf_le(&self, ray: &Ray, n_light: &Normal3f) -> Pdf {
        let it = Hit::new(
            ray.o,
            ray.time,
            Vector3f::ZERO,
            Vector3f::from(n_light),
            n_light.clone(),
            Some(self.medium_interface.clone()),
        );

        let pdf_pos = self.shape.pdf(&it);
        let pdf_dir = if self.two_sided {
            0.5 * cosine_hemisphere_pdf(n_light.abs_dot(&ray.d))
        } else {
            cosine_hemisphere_pdf(n_light.dot(&ray.d))
        };

        Pdf::new(pdf_pos, pdf_dir)
    }

    /// Returns the number of samples to use for the light source.
    fn get_num_samples(&self) -> usize {
        self.n_samples
    }

    /// Returns the area light's emitted radiance in a given outgoing direction.
    ///
    /// * `hit` - Point on a surface to evaluate emitted radiance.
    /// * `w`   - Outgoing direction.
    fn l(&self, hit: &Hit, w: &Vector3f) -> Spectrum {
        if self.two_sided || hit.n.dot(w) > 0.0 {
            self.l_emit
        } else {
            Spectrum::ZERO
        }
    }
}

impl From<(&ParamSet, ArcTransform, Option<ArcMedium>, ArcShape, usize)> for DiffuseAreaLight {
    /// Create a `DiffuseAreaLight` from given parameter set, light to world transform medium, shape and id.
    ///
    /// * `p` - A tuple containing the parameter set, light to world transform, medium, shape and id.
    fn from(p: (&ParamSet, ArcTransform, Option<ArcMedium>, ArcShape, usize)) -> Self {
        let (params, light_to_world, medium, shape, id) = p;

        let l = params.find_one_spectrum("L", Spectrum::ONE);
        let sc = params.find_one_spectrum("scale", Spectrum::ONE);
        let two_sided = params.find_one_bool("twosided", false);

        let mut n_samples = params.find_one_int("samples", params.find_one_int("nsamples", 1));
        if options().quick_render {
            n_samples = max(1, n_samples / 4);
        }

        Self::new(
            id,
            light_to_world,
            MediumInterface::from(medium),
            l * sc,
            n_samples as usize,
            shape,
            two_sided,
        )
    }
}
