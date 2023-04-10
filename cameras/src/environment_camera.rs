//! Environment Camera

use core::camera::*;
use core::film::*;
use core::geometry::*;
use core::medium::*;
use core::paramset::*;
use core::pbrt::*;
use std::mem::swap;
use std::sync::Arc;

// Environment camera.
pub struct EnvironmentCamera {
    /// Common camera parameters.
    pub data: CameraData,
}

impl EnvironmentCamera {
    /// Create a new environment camera.
    ///
    /// * `camera_to_world` - Animated transformation describing the camera's motion in the scene.
    /// * `shutter_open`    - Time when shutter is open.
    /// * `shutter_close`   - Time when shutter is closed.
    /// * `film`            - The film to capture the rendered image.
    /// * `medium`          - Scattering medium the camera lies in.
    pub fn new(
        camera_to_world: AnimatedTransform,
        shutter_open: Float,
        shutter_close: Float,
        film: Film,
        medium: Option<ArcMedium>,
    ) -> Self {
        Self {
            data: CameraData::new(camera_to_world, shutter_open, shutter_close, film, medium),
        }
    }
}

impl Camera for EnvironmentCamera {
    /// Returns the camera data.
    fn get_data(&self) -> &CameraData {
        &self.data
    }

    /// Returns a ray corresponding to a given sample. It also returns, a floating point value that affects how much the
    /// radiance arriving at the film plane will contribute to final image.
    ///
    /// * `sample` - The sample.
    fn generate_ray(&self, sample: &CameraSample) -> (Ray, Float) {
        // Compute environment camera ray direction.
        let theta = PI * sample.p_film.y / self.data.film.full_resolution.y as Float;
        let phi = TWO_PI * sample.p_film.x / self.data.film.full_resolution.x as Float;
        let dir = Vector3f::new(sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi));

        let ray = Ray::new(
            Point3f::new(0.0, 0.0, 0.0),
            dir,
            INFINITY,
            lerp(sample.time, self.data.shutter_open, self.data.shutter_close),
            self.data.medium.as_ref().map(Arc::clone),
        );

        (self.data.camera_to_world.transform_ray(&ray), 1.0)
    }

    /// Return the spatial and directional PDFs, as a tuple, for sampling a particular ray leaving the camera.
    ///
    /// * `ray` - The ray.
    fn pdf_we(&self, _ray: &Ray) -> PDFResult {
        panic!("NOT IMPLEMENTED");
    }
}

impl From<(&ParamSet, &AnimatedTransform, Film, Option<ArcMedium>)> for EnvironmentCamera {
    /// Create a `EnvironmentCamera` from given parameter set, animated transform, film and medium.
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

        Self::new(cam2world.clone(), shutter_open, shutter_close, film, medium)
    }
}
