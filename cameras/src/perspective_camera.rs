//! Perspective Camera

use core::camera::*;
use core::film::*;
use core::geometry::*;
use core::interaction::Interaction;
use core::interaction::SurfaceInteraction;
use core::light::VisibilityTester;
use core::medium::*;
use core::paramset::*;
use core::pbrt::*;
use core::sampling::*;
use core::spectrum::*;
use std::mem::swap;
use std::sync::Arc;

/// Perspective camera.
pub struct PerspectiveCamera {
    /// Common camera parameters.
    pub data: CameraData,

    /// Projective camera parameters.
    pub proj_data: ProjectiveCameraData,

    /// Differential change in x-coordinate of origin for camera rays.
    pub dx_camera: Vector3f,

    /// Differential change in y-coordinate of origin for camera rays.
    pub dy_camera: Vector3f,

    /// Area covered by the image plane bounds at z=1.
    pub a: Float,
}

impl PerspectiveCamera {
    /// Create a new perspective camera.
    ///
    /// * `camera_to_world` - Animated transformation describing the camera's
    ///                       motion in the scene.
    /// * `screen_window`   - Bounds of screen space.
    /// * `shutter_open`    - Time when shutter is open.
    /// * `shutter_close`   - Time when shutter is closed.
    /// * `lens_radius`     - Radius of camera lens.
    /// * `focal_distance`  - Focal distance.
    /// * `fov`             - The field-of-view angle in degrees.
    /// * `film`            - The film to capture the rendered image.
    /// * `medium`          - Scattering medium the camera lies in.
    pub fn new(
        camera_to_world: AnimatedTransform,
        screen_window: Bounds2f,
        shutter_open: Float,
        shutter_close: Float,
        lens_radius: Float,
        focal_distance: Float,
        fov: Float,
        film: Film,
        medium: Option<ArcMedium>,
    ) -> Self {
        let film_clone = film;
        let res = film_clone.full_resolution;

        let data = CameraData::new(
            camera_to_world,
            shutter_open,
            shutter_close,
            film_clone,
            medium,
        );
        let proj_data = ProjectiveCameraData::new(
            &data,
            Transform::perspective(fov, 1e-2, 1000.0),
            screen_window,
            lens_radius,
            focal_distance,
        );

        // Compute differential changes in origin for perspective camera rays.
        let dx_camera = proj_data
            .raster_to_camera
            .transform_point(&Point3f::new(1.0, 0.0, 0.0))
            - proj_data
                .raster_to_camera
                .transform_point(&Point3f::new(0.0, 0.0, 0.0));
        let dy_camera = proj_data
            .raster_to_camera
            .transform_point(&Point3f::new(0.0, 1.0, 0.0))
            - proj_data
                .raster_to_camera
                .transform_point(&Point3f::new(0.0, 0.0, 0.0));

        // Compute the image plane bounds at z=1 for perspective camera.
        let mut p_min = proj_data
            .raster_to_camera
            .transform_point(&Point3f::new(0.0, 0.0, 0.0));

        let mut p_max = proj_data.raster_to_camera.transform_point(&Point3f::new(
            res.x as Float,
            res.y as Float,
            0.0,
        ));

        p_min /= p_min.z;
        p_max /= p_max.z;

        let a = abs((p_max.x - p_min.x) * (p_max.y - p_min.y));

        Self {
            data,
            proj_data,
            dx_camera,
            dy_camera,
            a,
        }
    }
}

impl Camera for PerspectiveCamera {
    /// Returns the camera data.
    fn get_data(&self) -> &CameraData {
        &self.data
    }

    /// Returns a ray corresponding to a given sample. It also returns, a floating
    /// point value that affects how much the radiance arriving at the film plane
    /// will contribute to final image.
    ///
    /// * `sample` - The sample.
    fn generate_ray(&self, sample: &CameraSample) -> (Ray, Float) {
        // Compute raster and camera sample positions.
        let p_film = Point3f::new(sample.p_film.x, sample.p_film.y, 0.0);

        let p_camera = self.proj_data.raster_to_camera.transform_point(&p_film);

        let mut ray = Ray::new(
            Point3f::new(0.0, 0.0, 0.0),
            Vector3f::from(p_camera).normalize(),
            INFINITY,
            lerp(sample.time, self.data.shutter_open, self.data.shutter_close),
            self.data.medium.as_ref().map(Arc::clone),
        );

        // Modify ray for depth of field.
        if self.proj_data.lens_radius > 0.0 {
            // Sample point on lens.
            let p_lens = self.proj_data.lens_radius * concentric_sample_disk(&sample.p_lens);

            // Compute point on plane of focus.
            let ft = self.proj_data.focal_distance / ray.d.z;
            let p_focus = ray.at(ft);

            // Update ray for effect of lens.
            ray.o = Point3f::new(p_lens.x, p_lens.y, 0.0);
            ray.d = (p_focus - ray.o).normalize();
        }

        (self.data.camera_to_world.transform_ray(&ray), 1.0)
    }

    /// Returns a main ray and rays shifted one pixel in x and y directions on
    /// the film plane for corresponding to a given sample. It also returns a
    /// floating point value that affects how much the radiance arriving at the
    /// film plane will contribute to final image.
    ///
    /// * `sample` - The sample.
    fn generate_ray_differential(&self, sample: &CameraSample) -> (Ray, Float) {
        // Compute main perspective viewing ray.

        // Compute raster and camera sample positions.
        let p_film = Point3f::new(sample.p_film.x, sample.p_film.y, 0.0);

        let p_camera = self.proj_data.raster_to_camera.transform_point(&p_film);

        let mut ray = Ray::new(
            Point3f::new(0.0, 0.0, 0.0),
            Vector3f::from(p_camera).normalize(),
            INFINITY,
            lerp(sample.time, self.data.shutter_open, self.data.shutter_close),
            self.data.medium.as_ref().map(Arc::clone),
        );

        // Modify ray for depth of field.
        if self.proj_data.lens_radius > 0.0 {
            // Sample point on lens.
            let p_lens = self.proj_data.lens_radius * concentric_sample_disk(&sample.p_lens);

            // Compute point on plane of focus.
            let ft = self.proj_data.focal_distance / ray.d.z;
            let p_focus = ray.at(ft);

            // Update ray for effect of lens.
            ray.o = Point3f::new(p_lens.x, p_lens.y, 0.0);
            ray.d = (p_focus - ray.o).normalize();
        }

        // Compute ray differentials for perspective camera.
        let rd = if self.proj_data.lens_radius > 0.0 {
            // Compute perspective camera camera ray differentials accounting for lens.

            // Sample point on lens.
            let p_lens = self.proj_data.lens_radius * concentric_sample_disk(&sample.p_lens);

            let dx = Vector3f::from(p_camera + self.dx_camera).normalize();
            let ft = self.proj_data.focal_distance / dx.z;
            let p_focus = Point3f::new(0.0, 0.0, 0.0) + (ft * dx);
            let rx_origin = Point3f::new(p_lens.x, p_lens.y, 0.0);
            let rx_direction = (p_focus - rx_origin).normalize();

            let dy = Vector3f::from(p_camera + self.dy_camera).normalize();
            let ft = self.proj_data.focal_distance / dy.z;
            let p_focus = Point3f::new(0.0, 0.0, 0.0) + (ft * dy);
            let ry_origin = Point3f::new(p_lens.x, p_lens.y, 0.0);
            let ry_direction = (p_focus - ry_origin).normalize();

            RayDifferential::new(rx_origin, ry_origin, rx_direction, ry_direction)
        } else {
            let rx_origin = ray.o;
            let ry_origin = ray.o;
            let rx_direction = (Vector3f::from(p_camera) + self.dx_camera).normalize();
            let ry_direction = (Vector3f::from(p_camera) + self.dy_camera).normalize();
            RayDifferential::new(rx_origin, ry_origin, rx_direction, ry_direction)
        };
        ray.differentials = Some(rd);

        (self.data.camera_to_world.transform_ray(&ray), 1.0)
    }

    /// Evaluate the importance emitted from the point on the camera in a
    /// direction. The `include_raster_point` is true, then a raster position
    /// associated with the ray on the film is returned as well.
    ///
    /// * `ray`                  - The ray.
    /// * `include_raster_point` - Indicates whether or not to return the raster
    ///                            position.
    fn we(&self, ray: &Ray, include_raster: bool) -> (Spectrum, Option<Point2f>) {
        let mut p_raster2: Option<Point2f> = None;

        // Interpolate camera matrix and check if `ω` is forward-facing.
        let c2w = self.data.camera_to_world.interpolate(ray.time);
        let cos_theta = ray
            .d
            .dot(&c2w.transform_vector(&Vector3f::new(0.0, 0.0, 1.0)));
        if cos_theta <= 0.0 {
            return (Spectrum::ZERO, p_raster2);
        }

        // Map ray (p, ω) onto the raster grid.
        let p_focus = ray.at(if self.proj_data.lens_radius > 0.0 {
            self.proj_data.focal_distance
        } else {
            1.0
        } / cos_theta);

        let inv_raster_to_cam = self.proj_data.raster_to_camera.inverse();
        let inv_c2w = c2w.inverse();
        let p_raster = inv_raster_to_cam.transform_point(&inv_c2w.transform_point(&p_focus));

        // Return raster position if requested.
        if include_raster {
            p_raster2 = Some(Point2f::new(p_raster.x, p_raster.y));
        }

        // Return zero importance for out of bounds points.
        let sample_bounds = self.data.film.get_sample_bounds();
        if p_raster.x < sample_bounds.p_min.x as Float
            || p_raster.x >= sample_bounds.p_max.x as Float
            || p_raster.y < sample_bounds.p_min.y as Float
            || p_raster.y >= sample_bounds.p_max.y as Float
        {
            return (Spectrum::ZERO, p_raster2);
        }

        // Compute lens area of perspective camera.
        let lens_area = if self.proj_data.lens_radius != 0.0 {
            PI * self.proj_data.lens_radius * self.proj_data.lens_radius
        } else {
            1.0
        };

        // Return importance for point on image plane.
        let cos_2_theta = cos_theta * cos_theta;
        let importance = Spectrum::new(1.0 / (self.a * lens_area * cos_2_theta * cos_2_theta));
        (importance, p_raster2)
    }

    /// Return the spatial and directional PDFs, as a tuple, for sampling a
    /// particular ray leaving the camera.
    ///
    /// * `ray` - The ray.
    fn pdf_we(&self, ray: &Ray) -> PDFResult {
        // Interpolate camera matrix and check if `ω` is forward-facing.
        let c2w = self.data.camera_to_world.interpolate(ray.time);
        let cos_theta = ray
            .d
            .dot(&c2w.transform_vector(&Vector3f::new(0.0, 0.0, 1.0)));
        if cos_theta <= 0.0 {
            return PDFResult::default();
        }

        // Map ray (p, ω) onto the raster grid.
        let p_focus = ray.at(if self.proj_data.lens_radius > 0.0 {
            self.proj_data.focal_distance
        } else {
            1.0
        } / cos_theta);

        let inv_raster_to_cam = self.proj_data.raster_to_camera.inverse();
        let inv_c2w = c2w.inverse();
        let p_raster = inv_raster_to_cam.transform_point(&inv_c2w.transform_point(&p_focus));

        // Return zero probability for out of bounds points.
        let sample_bounds = self.data.film.get_sample_bounds();
        if p_raster.x < sample_bounds.p_min.x as Float
            || p_raster.x >= sample_bounds.p_max.x as Float
            || p_raster.y < sample_bounds.p_min.y as Float
            || p_raster.y >= sample_bounds.p_max.y as Float
        {
            return PDFResult::default();
        }

        // Compute lens area of perspective camera.
        let lens_area = if self.proj_data.lens_radius != 0.0 {
            PI * self.proj_data.lens_radius * self.proj_data.lens_radius
        } else {
            1.0
        };
        PDFResult::new(
            1.0 / lens_area,
            1.0 / (self.a * cos_theta * cos_theta * cos_theta),
        )
    }

    /// Returns a PDF value with respect to the solid angle at a reference point.
    ///
    /// * `interaction`          - The interaction point.
    /// * `u`                    - Used to sample point on the lens.
    /// * `include_raster_point` - Indicates whether or not to return the raster
    ///                            position.
    fn sample_wi(&self, interaction: &Interaction, u: &Point2f) -> SampleResult {
        // Uniformly sample a lens interaction `lensIntr`.
        let interaction_hit = interaction.get_hit().to_owned();
        let time = interaction_hit.time;

        let p_lens = self.proj_data.lens_radius * concentric_sample_disk(u);
        let p_lens_world = self
            .data
            .camera_to_world
            .transform_point(time, &Point3f::new(p_lens.x, p_lens.y, 0.0));

        let mut si = SurfaceInteraction::new(
            p_lens_world,
            Vector3f::default(),
            Point2f::default(),
            Vector3f::default(),
            Vector3f::default(),
            Vector3f::default(),
            Normal3f::default(),
            Normal3f::default(),
            time,
            None,
            0,
        );
        si.hit.n = Normal3f::from(
            &self
                .data
                .camera_to_world
                .transform_vector(time, &Vector3f::new(0.0, 0.0, 1.0)),
        );
        let lens_intr = Interaction::Surface { si };
        let lens_intr_hit = lens_intr.get_hit().to_owned();

        let lens_intr_hit_p = lens_intr_hit.p;
        let lens_intr_hit_n = lens_intr_hit.n;
        let interaction_hit_p = interaction_hit.p;

        // Populate arguments and compute the importance value.
        let vis = VisibilityTester::new(interaction_hit, lens_intr_hit);
        let mut wi = lens_intr_hit_p - interaction_hit_p;
        let dist = wi.length();
        wi /= dist;

        // Compute PDF for importance arriving at `interaction`.

        // Compute lens area of perspective camera.
        let lens_area = if self.proj_data.lens_radius != 0.0 {
            PI * self.proj_data.lens_radius * self.proj_data.lens_radius
        } else {
            1.0
        };
        let pdf = (dist * dist) / (lens_intr_hit_n.abs_dot(&wi) * lens_area);
        let (spectrum, p_raster) = self.we(&lens_intr.spawn_ray(&(-wi)), true);

        SampleResult::new(spectrum, wi, pdf, p_raster, vis)
    }
}

impl From<(&ParamSet, &AnimatedTransform, Film, Option<ArcMedium>)> for PerspectiveCamera {
    /// Create a `PerspectiveCamera` from given parameter set, animated transform,
    /// film and medium.
    ///
    /// * `p` - A tuple containing  parameter set, animated transform, film and
    ///         medium.
    fn from(p: (&ParamSet, &AnimatedTransform, Film, Option<ArcMedium>)) -> Self {
        let (params, cam2world, film, medium) = p;

        // Extract common camera parameters from `ParamSet`
        let mut shutter_open = params.find_one_float("shutteropen", 0.0);
        let mut shutter_close = params.find_one_float("shutterclose", 1.0);
        if shutter_close < shutter_open {
            warn!(
                "Shutter close time [{}] < shutter open [{}]. 
                Swapping them.",
                shutter_close, shutter_open
            );
            swap(&mut shutter_close, &mut shutter_open);
        }

        let lens_radius = params.find_one_float("lensradius", 0.0);
        let focal_distance = params.find_one_float("focaldistance", 1e6);

        let frame = params.find_one_float(
            "frameaspectratio",
            film.full_resolution.x as Float / film.full_resolution.y as Float,
        );
        let mut screen = if frame > 1.0 {
            Bounds2::new(Point2::new(-frame, -1.0), Point2::new(frame, 1.0))
        } else {
            Bounds2::new(
                Point2::new(-1.0, -1.0 / frame),
                Point2::new(1.0, 1.0 / frame),
            )
        };

        let sw = params.find_float("screenwindow");
        let swi = sw.len();
        if swi > 0 {
            if swi == 4 {
                screen.p_min.x = sw[0];
                screen.p_max.x = sw[1];
                screen.p_min.y = sw[2];
                screen.p_max.y = sw[3];
            } else {
                error!("'screenwindow' should have four values");
            }
        }

        let mut fov = params.find_one_float("fov", 90.0);
        let half_fov = params.find_one_float("halffov", -1.0);
        if half_fov > 0.0 {
            // Hack for structure synth, which exports half of the full
            // fov.
            fov = 2.0 * half_fov;
        }

        Self::new(
            cam2world.clone(),
            screen,
            shutter_open,
            shutter_close,
            lens_radius,
            focal_distance,
            fov,
            film,
            medium,
        )
    }
}
