//! Orthographic Camera

use core::camera::*;
use core::film::*;
use core::geometry::*;
use core::medium::*;
use core::paramset::*;
use core::pbrt::*;
use core::sampling::*;
use std::mem::swap;
use std::sync::Arc;

/// Orthographic camera.
pub struct OrthographicCamera {
    /// Common camera parameters.
    pub data: CameraData,

    /// Projective camera parameters.
    pub proj_data: ProjectiveCameraData,

    /// Differential change in x-coordinate of origin for camera rays.
    pub dx_camera: Vector3f,

    /// Differential change in y-coordinate of origin for camera rays.
    pub dy_camera: Vector3f,
}

impl OrthographicCamera {
    /// Create a new orthographic camera.
    ///
    /// * `camera_to_world` - Animated transformation describing the camera's motion in the scene.
    /// * `screen_window`   - Bounds of screen space.
    /// * `shutter_open`    - Time when shutter is open.
    /// * `shutter_close`   - Time when shutter is closed.
    /// * `lens_radius`     - Radius of camera lens.
    /// * `focal_distance`  - Focal distance.
    /// * `film`            - The film to capture the rendered image.
    /// * `medium`          - Scattering medium the camera lies in.
    pub fn new(
        camera_to_world: AnimatedTransform,
        screen_window: Bounds2f,
        shutter_open: Float,
        shutter_close: Float,
        lens_radius: Float,
        focal_distance: Float,
        film: Film,
        medium: Option<ArcMedium>,
    ) -> Self {
        let data = CameraData::new(camera_to_world, shutter_open, shutter_close, film, medium);
        let proj_data = ProjectiveCameraData::new(
            &data,
            Transform::orthographic(0.0, 1.0),
            screen_window,
            lens_radius,
            focal_distance,
        );

        // Compute differential changes in origin for orthographic camera rays.
        let dx_camera = proj_data
            .raster_to_camera
            .transform_vector(&Vector3f::new(1.0, 0.0, 0.0));
        let dy_camera = proj_data
            .raster_to_camera
            .transform_vector(&Vector3f::new(0.0, 1.0, 0.0));

        Self {
            data,
            proj_data,
            dx_camera,
            dy_camera,
        }
    }
}

impl Camera for OrthographicCamera {
    /// Returns the camera data.
    fn get_data(&self) -> &CameraData {
        &self.data
    }

    /// Returns a ray corresponding to a given sample. It also returns, a floating point value that affects how much the
    /// radiance arriving at the film plane will contribute to final image.
    ///
    /// * `sample` - The sample.
    fn generate_ray(&self, sample: &CameraSample) -> (Ray, Float) {
        // Compute raster and camera sample positions.
        let p_film = Point3f::new(sample.p_film.x, sample.p_film.y, 0.0);

        let p_camera = self.proj_data.raster_to_camera.transform_point(&p_film);

        let mut ray = Ray::new(
            p_camera,
            Vector3f::new(0.0, 0.0, 1.0),
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

    /// Returns a main ray and rays shifted one pixel in x and y directions on the film plane for corresponding to a
    /// given sample. It also returns a floating point value that affects how much the radiance arriving at the film
    /// plane will contribute to final image.
    ///
    /// * `sample` - The sample.
    fn generate_ray_differential(&self, sample: &CameraSample) -> (Ray, Float) {
        // Compute main orthographic viewing ray.

        // Compute raster and camera sample positions.
        let p_film = Point3f::new(sample.p_film.x, sample.p_film.y, 0.0);

        let p_camera = self.proj_data.raster_to_camera.transform_point(&p_film);

        let mut ray = Ray::new(
            p_camera,
            Vector3f::new(0.0, 0.0, 1.0),
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

        // Compute ray differentials for orthographic camera.
        let rd = if self.proj_data.lens_radius > 0.0 {
            // Compute orthographic camera ray differentials accounting for lens.

            // Sample point on lens.
            let p_lens = self.proj_data.lens_radius * concentric_sample_disk(&sample.p_lens);
            let ft = self.proj_data.focal_distance / ray.d.z;

            let p_focus = p_camera + self.dx_camera + (ft * Vector3f::new(0.0, 0.0, 1.0));
            let rx_origin = Point3f::new(p_lens.x, p_lens.y, 0.0);
            let rx_direction = (p_focus - rx_origin).normalize();

            let p_focus = p_camera + self.dy_camera + (ft * Vector3f::new(0.0, 0.0, 1.0));
            let ry_origin = Point3f::new(p_lens.x, p_lens.y, 0.0);
            let ry_direction = (p_focus - ry_origin).normalize();

            RayDifferential::new(rx_origin, ry_origin, rx_direction, ry_direction)
        } else {
            let rx_origin = ray.o + self.dx_camera;
            let ry_origin = ray.o + self.dy_camera;
            let rx_direction = ray.d;
            let ry_direction = ray.d;
            RayDifferential::new(rx_origin, ry_origin, rx_direction, ry_direction)
        };
        ray.differentials = Some(rd);

        (self.data.camera_to_world.transform_ray(&ray), 1.0)
    }

    /// Return the spatial and directional PDFs, as a tuple, for sampling a particular ray leaving the camera.
    ///
    /// * `ray` - The ray.
    fn pdf_we(&self, _ray: &Ray) -> PDFResult {
        panic!("NOT IMPLEMENTED");
    }
}

impl From<(&ParamSet, &AnimatedTransform, Film, Option<ArcMedium>)> for OrthographicCamera {
    /// Create a `OrthographicCamera` from given parameter set, animated transform, film and medium.
    ///
    /// * `p` - A tuple containing  parameter set, animated transform, film and medium.
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
            Bounds2::new(Point2::new(-1.0, -1.0 / frame), Point2::new(1.0, 1.0 / frame))
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

        Self::new(
            cam2world.clone(),
            screen,
            shutter_open,
            shutter_close,
            lens_radius,
            focal_distance,
            film,
            medium,
        )
    }
}
