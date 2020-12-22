//! Orthographic Camera

#![allow(dead_code)]
use super::{
    concentric_sample_disk, lerp, AnimatedTransform, ArcMedium, Bounds2f, Camera, CameraData,
    CameraSample, Film, Float, PDFResult, Point3f, ProjectiveCameraData, Ray, RayDifferential,
    Transform, Vector3f, INFINITY,
};
use std::sync::Arc;

/// Orthographic camera.
#[derive(Clone)]
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
    /// * `camera_to_world` - Animated transformation describing the camera's
    ///                       motion in the scene.
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
        film: Arc<Film>,
        medium: ArcMedium,
    ) -> Self {
        let data = CameraData::new(
            camera_to_world,
            shutter_open,
            shutter_close,
            film.clone(),
            medium.clone(),
        );
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
            p_camera,
            Vector3f::new(0.0, 0.0, 1.0),
            INFINITY,
            lerp(sample.time, self.data.shutter_open, self.data.shutter_close),
            Some(self.data.medium.clone()),
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
        // Compute main orthographic viewing ray.

        // Compute raster and camera sample positions.
        let p_film = Point3f::new(sample.p_film.x, sample.p_film.y, 0.0);

        let p_camera = self.proj_data.raster_to_camera.transform_point(&p_film);

        let mut ray = Ray::new(
            p_camera,
            Vector3f::new(0.0, 0.0, 1.0),
            INFINITY,
            lerp(sample.time, self.data.shutter_open, self.data.shutter_close),
            Some(self.data.medium.clone()),
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

    /// Return the spatial and directional PDFs, as a tuple, for sampling a
    /// particular ray leaving the camera.
    ///
    /// * `ray` - The ray.
    fn pdf_we(&self, _ray: &Ray) -> PDFResult {
        panic!("NOT IMPLEMENTED");
    }
}
