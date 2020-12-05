//! Perspective Camera

#![allow(dead_code)]
use super::{
    abs, concentric_sample_disk, lerp, AnimatedTransform, ArcMedium, Bounds2f, Camera, CameraData,
    CameraSample, Film, Float, PDFResult, Point3f, ProjectiveCameraData, Ray, RayDifferential,
    Transform, Vector3f, INFINITY,
};

/// Perspective camera.
#[derive(Clone)]
pub struct PerspectiveCamera {
    /// Common camera parameters.
    pub camera_data: CameraData,

    /// Projective camera parameters.
    pub projective_data: ProjectiveCameraData,

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
        medium: ArcMedium,
    ) -> Self {
        let camera_data = CameraData::new(
            camera_to_world,
            shutter_open,
            shutter_close,
            film,
            medium.clone(),
        );
        let projective_data = ProjectiveCameraData::new(
            &camera_data,
            Transform::perspective(fov, 1e-2, 1000.0),
            screen_window,
            lens_radius,
            focal_distance,
        );

        // Compute differential changes in origin for perspective camera rays.
        let dx_camera = projective_data
            .raster_to_camera
            .transform_point(&Point3f::new(1.0, 0.0, 0.0))
            - projective_data
                .raster_to_camera
                .transform_point(&Point3f::new(0.0, 0.0, 0.0));
        let dy_camera = projective_data
            .raster_to_camera
            .transform_point(&Point3f::new(0.0, 1.0, 0.0))
            - projective_data
                .raster_to_camera
                .transform_point(&Point3f::new(0.0, 0.0, 0.0));

        // Compute the image plane bounds at z=1 for perspective camera.
        let res = film.full_resolution;

        let mut p_min = projective_data
            .raster_to_camera
            .transform_point(&Point3f::new(0.0, 0.0, 0.0));
        p_min /= p_min.z;

        let mut p_max = projective_data
            .raster_to_camera
            .transform_point(&Point3f::new(res.x as Float, res.y as Float, 0.0));
        p_max /= p_max.z;

        let a = abs((p_max.x - p_min.x) * (p_max.y - p_min.y));

        Self {
            camera_data,
            projective_data,
            dx_camera,
            dy_camera,
            a,
        }
    }
}

impl Camera for PerspectiveCamera {
    /// Returns a ray corresponding to a given sample. It also returns, a floating
    /// point value that affects how much the radiance arriving at the film plane
    /// will contribute to final image.
    ///
    /// * `sample` - The sample.
    fn generate_ray(&self, sample: &CameraSample) -> (Ray, Float) {
        // Compute raster and camera sample positions.
        let p_film = Point3f::new(sample.p_film.x, sample.p_film.y, 0.0);

        let p_camera = self
            .projective_data
            .raster_to_camera
            .transform_point(&p_film);

        let mut ray = Ray::new(
            Point3f::new(0.0, 0.0, 0.0),
            Vector3f::from(p_camera).normalize(),
            INFINITY,
            lerp(
                sample.time,
                self.camera_data.shutter_open,
                self.camera_data.shutter_close,
            ),
            Some(self.camera_data.medium.clone()),
        );

        // Modify ray for depth of field.
        if self.projective_data.lens_radius > 0.0 {
            // Sample point on lens.
            let p_lens = self.projective_data.lens_radius * concentric_sample_disk(&sample.p_lens);

            // Compute point on plane of focus.
            let ft = self.projective_data.focal_distance / ray.d.z;
            let p_focus = ray.at(ft);

            // Update ray for effect of lens.
            ray.o = Point3f::new(p_lens.x, p_lens.y, 0.0);
            ray.d = (p_focus - ray.o).normalize();
        }

        (self.camera_data.camera_to_world.transform_ray(&ray), 1.0)
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

        let p_camera = self
            .projective_data
            .raster_to_camera
            .transform_point(&p_film);

        let mut ray = Ray::new(
            Point3f::new(0.0, 0.0, 0.0),
            Vector3f::from(p_camera).normalize(),
            INFINITY,
            lerp(
                sample.time,
                self.camera_data.shutter_open,
                self.camera_data.shutter_close,
            ),
            Some(self.camera_data.medium.clone()),
        );

        // Modify ray for depth of field.
        if self.projective_data.lens_radius > 0.0 {
            // Sample point on lens.
            let p_lens = self.projective_data.lens_radius * concentric_sample_disk(&sample.p_lens);

            // Compute point on plane of focus.
            let ft = self.projective_data.focal_distance / ray.d.z;
            let p_focus = ray.at(ft);

            // Update ray for effect of lens.
            ray.o = Point3f::new(p_lens.x, p_lens.y, 0.0);
            ray.d = (p_focus - ray.o).normalize();
        }

        // Compute ray differentials for perspective camera.
        let rd = if self.projective_data.lens_radius > 0.0 {
            // Compute perspective camera camera ray differentials accounting for lens.

            // Sample point on lens.
            let p_lens = self.projective_data.lens_radius * concentric_sample_disk(&sample.p_lens);

            let dx = Vector3f::from(p_camera + self.dx_camera).normalize();
            let ft = self.projective_data.focal_distance / dx.z;
            let p_focus = Point3f::new(0.0, 0.0, 0.0) + (ft * dx);
            let rx_origin = Point3f::new(p_lens.x, p_lens.y, 0.0);
            let rx_direction = (p_focus - rx_origin).normalize();

            let dy = Vector3f::from(p_camera + self.dy_camera).normalize();
            let ft = self.projective_data.focal_distance / dy.z;
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

        (self.camera_data.camera_to_world.transform_ray(&ray), 1.0)
    }

    /// Return the spatial and directional PDFs, as a tuple, for sampling a
    /// particular ray leaving the camera.
    ///
    /// * `ray` - The ray.
    fn pdf_we(&self, _ray: &Ray) -> PDFResult {
        panic!("NOT IMPLEMENTED");
    }
}
