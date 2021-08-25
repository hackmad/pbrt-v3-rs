//! Whitted Integrator

#![allow(dead_code)]

use core::camera::*;
use core::geometry::*;
use core::integrator::*;
use core::light::*;
use core::material::*;
use core::paramset::*;
use core::pbrt::*;
use core::reflection::*;
use core::sampler::*;
use core::scene::*;
use core::spectrum::*;
use itertools::iproduct;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};

/// Implements Whitted's ray tracing algorithm.
pub struct WhittedIntegrator {
    /// Sampler responsible for choosing points on the image plane from which
    /// to trace rays and for supplying sample positions used by integrators.
    sampler: ArcSampler,

    /// The camera.
    camera: Arc<Mutex<ArcCamera>>,

    /// Pixel bounds for the image.
    pixel_bounds: Bounds2i,

    /// Maximum recursion depth.
    max_depth: usize,
}

impl WhittedIntegrator {
    /// Create a new `WhittedIntegrator`.
    ///
    /// * `max_depth`    - Maximum recursion depth.
    /// * `camera`       - The camera.
    /// * `sampler`      - The sampler.
    /// * `pixel_bounds` - Pixel bounds for the image.
    pub fn new(
        max_depth: usize,
        camera: ArcCamera,
        sampler: ArcSampler,
        pixel_bounds: Bounds2i,
    ) -> Self {
        Self {
            camera: Arc::new(Mutex::new(Arc::clone(&camera))),
            sampler,
            pixel_bounds,
            max_depth,
        }
    }

    /// Trace rays for specular reflection.
    ///
    /// * `ray`     - The ray.
    /// * `isect`   - The surface interaction.
    /// * `scene`   - The scene.
    /// * `sampler` - The sampler.
    /// * `depth`   - The recursive depth.
    fn specular_reflect(
        &self,
        ray: &mut Ray,
        isect: &SurfaceInteraction,
        scene: Arc<Scene>,
        sampler: &mut ArcSampler,
        depth: usize,
    ) -> Spectrum {
        if let Some(bsdf) = isect.bsdf.clone() {
            // Compute specular reflection direction `wi` and BSDF value.
            let wo = isect.hit.wo;

            let samp = Arc::get_mut(sampler).unwrap();
            let sample = samp.get_2d();
            let bxdf_type = BxDFType::from(BSDF_REFLECTION | BSDF_SPECULAR);
            let BxDFSample {
                f,
                pdf,
                wi,
                sampled_type: _,
            } = bsdf.sample_f(&wo, &sample, bxdf_type);

            // Return contribution of specular reflection
            let ns = isect.shading.n;
            if pdf > 0.0 && !f.is_black() && wi.abs_dot(&ns) != 0.0 {
                // Compute ray differential `rd` for specular reflection.
                let mut rd = isect.hit.spawn_ray(&wi);
                if let Some(differentials) = ray.differentials {
                    let rx_origin = isect.hit.p + isect.dpdx;
                    let ry_origin = isect.hit.p + isect.dpdy;

                    // Compute differential reflected directions.
                    let dndx = isect.shading.dndu * isect.dudx + isect.shading.dndv * isect.dvdx;
                    let dndy = isect.shading.dndu * isect.dudy + isect.shading.dndv * isect.dvdy;
                    let dwodx = -differentials.rx_direction - wo;
                    let dwody = -differentials.ry_direction - wo;
                    let ddndx = dwodx.dot(&ns) + wo.dot(&dndx);
                    let ddndy = dwody.dot(&ns) + wo.dot(&dndy);
                    let rx_direction =
                        wi - dwodx + 2.0 * Vector3f::from(wo.dot(&ns) * dndx + ddndx * ns);
                    let ry_direction =
                        wi - dwody + 2.0 * Vector3f::from(wo.dot(&ns) * dndy + ddndy * ns);
                    rd.differentials = Some(RayDifferential::new(
                        rx_origin,
                        ry_origin,
                        rx_direction,
                        ry_direction,
                    ));
                }

                return f
                    * self.li(&mut rd, Arc::clone(&scene), sampler, depth + 1)
                    * wi.abs_dot(&ns)
                    / pdf;
            }
        }

        Spectrum::new(0.0)
    }

    /// Trace rays for specular refraction.
    ///
    /// * `ray`     - The ray.
    /// * `isect`   - The surface interaction.
    /// * `scene`   - The scene.
    /// * `sampler` - The sampler.
    /// * `depth`   - The recursive depth.
    fn specular_transmit(
        &self,
        ray: &mut Ray,
        isect: &SurfaceInteraction,
        scene: Arc<Scene>,
        sampler: &mut ArcSampler,
        depth: usize,
    ) -> Spectrum {
        if let Some(bsdf) = isect.bsdf.clone() {
            let wo = isect.hit.wo;
            let p = isect.hit.p;

            let samp = Arc::get_mut(sampler).unwrap();
            let sample = samp.get_2d();
            let bxdf_type = BxDFType::from(BSDF_TRANSMISSION | BSDF_SPECULAR);
            let BxDFSample {
                f,
                pdf,
                wi,
                sampled_type: _,
            } = bsdf.sample_f(&wo, &sample, bxdf_type);

            let mut ns = isect.shading.n;
            if pdf > 0.0 && !f.is_black() && wi.abs_dot(&ns) != 0.0 {
                // Compute ray differential _rd_ for specular transmission
                let mut rd = isect.hit.spawn_ray(&wi);
                if let Some(differentials) = ray.differentials {
                    let rx_origin = p + isect.dpdx;
                    let ry_origin = p + isect.dpdy;

                    let mut dndx =
                        isect.shading.dndu * isect.dudx + isect.shading.dndv * isect.dvdx;
                    let mut dndy =
                        isect.shading.dndu * isect.dudy + isect.shading.dndv * isect.dvdy;

                    // The BSDF stores the IOR of the interior of the object being
                    // intersected. Compute the relative IOR by first out by assuming
                    // that the ray is entering the object.
                    let mut eta = 1.0 / bsdf.eta;
                    if wo.dot(&ns) < 0.0 {
                        // If the ray isn't entering, then we need to invert the
                        // relative IOR and negate the normal and its derivatives.
                        eta = 1.0 / eta;
                        ns = -ns;
                        dndx = -dndx;
                        dndy = -dndy;
                    }

                    /*
                      Notes on the derivation:
                      - pbrt computes the refracted ray as:
                        wi = -eta * omega_o + [ eta * (wo dot N) - cos(theta_t) ] * N
                        It flips the normal to lie in the same hemisphere as wo,
                        and then eta is the relative IOR from wo's medium to wi's
                        medium.

                      - If we denote the term in brackets by mu, then we have:
                        wi = -eta * omega_o + mu * N

                      - Now let's take the partial derivative.
                        We get: -eta * d/dx(omega_o) + mu * d/dx(N) + d/dx(mu) * N.

                      - We have the values of all of these except for d/dx(mu)
                        (using bits from the derivation of specularly reflected
                        ray deifferentials).

                      - The first term of d/dx(mu) is easy: eta d/dx(wo dot N).
                        We already have d/dx(wo dot N).

                      - The second term takes a little more work. We have:
                        cos(theta_i) = sqrt[1 - eta^2 * (1 - (wo dot N)^2)].

                        Starting from (wo dot N)^2 and reading outward,
                        we have cos^2(theta_o),
                        then sin^2(theta_o),
                        then sin^2(theta_i) (via Snell's law),
                        then cos^2(theta_i) and then cos(theta_i).

                      - Let's take the partial derivative of the sqrt expression.
                        We get:
                        (1 / 2) * (1 / cos(theta_i) * d/dx(1 - eta^2 * (1 - (wo dot N)^2)))

                      - That partial derivatve is equal to:
                        d/dx(eta^2 * (wo dot N)^2) = 2 * eta^2 * (wo dot N) * d/dx(wo dot N)

                      - Plugging it in, we have d/dx(mu) =
                        eta * d/dx(wo dot N) - (eta^2 * (wo dot N) * d/dx(wo dot N))/(-wi dot N)
                    */
                    let dwodx = -differentials.rx_direction - wo;
                    let dwody = -differentials.ry_direction - wo;
                    let ddndx = dwodx.dot(&ns) + wo.dot(&dndx);
                    let ddndy = dwody.dot(&ns) + wo.dot(&dndy);

                    let mu = eta * wo.dot(&ns) - wi.abs_dot(&ns);
                    let dmudx = (eta - (eta * eta * wo.dot(&ns)) / wi.abs_dot(&ns)) * ddndx;
                    let dmudy = (eta - (eta * eta * wo.dot(&ns)) / wi.abs_dot(&ns)) * ddndy;

                    let rx_direction = wi - eta * dwodx + Vector3f::from(mu * dndx + dmudx * ns);
                    let ry_direction = wi - eta * dwody + Vector3f::from(mu * dndy + dmudy * ns);

                    rd.differentials = Some(RayDifferential::new(
                        rx_origin,
                        ry_origin,
                        rx_direction,
                        ry_direction,
                    ));
                }

                return f
                    * self.li(&mut rd, Arc::clone(&scene), sampler, depth + 1)
                    * wi.abs_dot(&ns)
                    / pdf;
            }
        }

        Spectrum::new(0.0)
    }
}

impl Integrator for WhittedIntegrator {
    /// Render the scene.
    ///
    /// * `scene`        - The scene.
    fn render(&mut self, scene: Arc<Scene>) {
        // Compute number of tiles, `n_tiles`, to use for parallel rendering.
        let sample_bounds = Arc::clone(&self.camera)
            .lock()
            .unwrap()
            .get_film_sample_bounds();
        let sample_extent = sample_bounds.diagonal();
        let tile_size = 16;
        let n_tiles = Point2::new(
            ((sample_extent.x + tile_size - 1) / tile_size) as usize,
            ((sample_extent.y + tile_size - 1) / tile_size) as usize,
        );

        info!("Rendering {}x{} tiles", n_tiles.x, n_tiles.y);

        // Parallelize.
        let tiles = iproduct!(0..n_tiles.x, 0..n_tiles.y).par_bridge();
        tiles.for_each(|(tile_x, tile_y)| {
            let camera_clone = Arc::clone(&self.camera);

            // Render section of image corresponding to `tile`.
            let tile = Point2::new(tile_x, tile_y);

            // Get sampler instance for tile.
            let seed = tile.y * n_tiles.x + tile.x;
            let mut tile_sampler = Sampler::clone(&*self.sampler, seed as u64);

            let samples_per_pixel = {
                let tile_sampler_data = Arc::get_mut(&mut tile_sampler).unwrap().get_data();
                tile_sampler_data.samples_per_pixel
            };

            // Compute sample bounds for tile.
            let x0 = sample_bounds.p_min.x + tile.x as i32 * tile_size;
            let x1 = min(x0 + tile_size, sample_bounds.p_max.x);
            let y0 = sample_bounds.p_min.y + tile.y as i32 * tile_size;
            let y1 = min(y0 + tile_size, sample_bounds.p_max.y);
            let tile_bounds = Bounds2i::new(Point2i::new(x0, y0), Point2i::new(x1, y1));

            info!(
                "Starting image tile ({}, {}) -> {:}",
                tile_x, tile_y, tile_bounds
            );

            // Get `FilmTile` for tile.
            let mut film_tile = {
                let camera = camera_clone.lock().unwrap();
                camera.get_film_tile(tile_bounds)
            };

            // Loop over pixels in tile to render them.
            for pixel in tile_bounds {
                Arc::get_mut(&mut tile_sampler).unwrap().start_pixel(&pixel);

                // Do this check after the StartPixel() call; this keeps the
                // usage of RNG values from (most) Samplers that use RNGs
                // consistent, which improves reproducability / debugging.
                if !self.pixel_bounds.contains_exclusive(&pixel) {
                    continue;
                }

                loop {
                    // Initialize `CameraSample` for current sample.
                    let camera_sample = Arc::get_mut(&mut tile_sampler)
                        .unwrap()
                        .get_camera_sample(&pixel);

                    // Generate camera ray for current sample.
                    let (mut ray, ray_weight) = {
                        let camera = camera_clone.lock().unwrap();
                        camera.generate_ray_differential(&camera_sample)
                    };
                    ray.scale_differentials(1.0 / (samples_per_pixel as Float).sqrt());

                    // Evaluate radiance along camera ray.
                    let mut l = Spectrum::new(0.0);
                    if ray_weight > 0.0 {
                        l = self.li(&mut ray, scene.clone(), &mut tile_sampler, 0);
                    }

                    // Issue warning if unexpected radiance value returned.
                    let tile_sampler_data = Arc::get_mut(&mut tile_sampler).unwrap().get_data();
                    let current_sample_number = tile_sampler_data.current_sample_number();
                    if l.has_nans() {
                        error!(
                            "Not-a-number radiance value returned for pixel
                                ({}, {}), sample {}. Setting to black.",
                            pixel.x, pixel.y, current_sample_number
                        );
                        l = Spectrum::new(0.0);
                    } else if l.y() < -1e-5 {
                        error!(
                            "Negative luminance value, {}, returned for pixel
                                ({}, {}), sample {}. Setting to black.",
                            l.y(),
                            pixel.x,
                            pixel.y,
                            current_sample_number
                        );
                        l = Spectrum::new(0.0);
                    } else if l.y().is_infinite() {
                        error!(
                            "Infinite luminance value returned for pixel
                                ({}, {}), sample {}. Setting to black.",
                            pixel.x, pixel.y, current_sample_number
                        );
                        l = Spectrum::new(0.0);
                    }

                    debug!(
                        "Pixel: {:}, Camera sample: {:} -> ray: {:}, ray weight {} -> L = {:}",
                        pixel, camera_sample, ray, ray_weight, l
                    );

                    // Add camera ray's contribution to image.
                    film_tile.add_sample(camera_sample.p_film, l, ray_weight);

                    if !Arc::get_mut(&mut tile_sampler).unwrap().start_next_sample() {
                        break;
                    }
                }
            }

            info!(
                "Finished image tile ({}, {}) -> {:}",
                tile_x, tile_y, tile_bounds
            );

            // Merge image tile into `Film`.
            let mut camera = camera_clone.lock().unwrap();
            Arc::get_mut(&mut *camera)
                .unwrap()
                .merge_film_tile(&film_tile);
        });

        info!("Rendering finished.");

        // Save final image after rendering.
        let camera_clone = Arc::clone(&self.camera);
        let mut camera = camera_clone.lock().unwrap();
        Arc::get_mut(&mut *camera).unwrap().write_image(1.0);
        info!("Output image written.");
    }

    /// Returns the incident radiance at the origin of a given ray.
    ///
    /// * `ray`     - The ray.
    /// * `scene`   - The scene.
    /// * `sampler` - The sampler.
    /// * `depth`   - The recursion depth.
    fn li(
        &self,
        ray: &mut Ray,
        scene: Arc<Scene>,
        sampler: &mut ArcSampler,
        depth: usize,
    ) -> Spectrum {
        let mut l = Spectrum::new(0.0);

        // Find closest ray intersection or return background radiance.
        if let Some(mut isect) = scene.intersect(ray) {
            // Compute emitted and reflected light at ray intersection point.

            // Initialize common variables for Whitted integrator.
            let n = isect.shading.n;
            let wo = isect.hit.wo;

            // Compute scattering functions for surface interaction.
            isect.compute_scattering_functions(ray, false, TransportMode::Radiance);
            if isect.bsdf.is_none() {
                let mut new_ray = isect.hit.spawn_ray(&ray.d);
                return self.li(&mut new_ray, scene.clone(), sampler, depth);
            }

            // Compute emitted light if ray hit an area light source.
            l += isect.le(&wo);

            // Add contribution of each light source.
            for light in scene.lights.iter() {
                let sample = Arc::get_mut(sampler).unwrap().get_2d();
                let Li {
                    wi,
                    pdf,
                    visibility,
                    value: li,
                } = light.sample_li(&isect.hit, &sample);

                if li.is_black() || pdf == 0.0 {
                    continue;
                }

                let f = isect
                    .bsdf
                    .clone()
                    .unwrap()
                    .f(&wo, &wi, BxDFType::from(BSDF_ALL));

                // If no visiblity tester, then unoccluded = true.
                let unoccluded = visibility.map_or(true, |vis| vis.unoccluded(scene.clone()));
                if !f.is_black() && unoccluded {
                    l += f * li * wi.abs_dot(&n) / pdf;
                }
            }
            if depth + 1 < self.max_depth {
                // Trace rays for specular reflection and refraction.
                l += self.specular_reflect(ray, &isect, Arc::clone(&scene), sampler, depth);
                l += self.specular_transmit(ray, &isect, Arc::clone(&scene), sampler, depth);
            }
        } else {
            if let Some(rd) = ray.differentials {
                for light in scene.lights.iter() {
                    l += light.le(&rd);
                }
            }
        }

        l
    }
}

impl From<(&ParamSet, ArcSampler, ArcCamera)> for WhittedIntegrator {
    /// Create a `WhittedIntegrator` from given parameter set and camera.
    ///
    /// * `p` - A tuple containing parameter set and camera.
    fn from(p: (&ParamSet, ArcSampler, ArcCamera)) -> Self {
        let (params, sampler, camera) = p;

        let max_depth = params.find_one_int("max_depth", 5) as usize;

        let pb = params.find_int("pixelbounds");
        let np = pb.len();

        let mut pixel_bounds = camera.get_film_sample_bounds();
        if np > 0 {
            if np != 4 {
                error!("Expected 4 values for 'pixel_bounds' parameter. Got {}", np);
            } else {
                pixel_bounds = pixel_bounds.intersect(&Bounds2i::new(
                    Point2i::new(pb[0], pb[1]),
                    Point2i::new(pb[2], pb[3]),
                ));
                if pixel_bounds.area() == 0 {
                    error!("Degenerate 'pixel_bounds' specified.");
                }
            }
        }

        Self::new(
            max_depth,
            Arc::clone(&camera),
            Arc::clone(&sampler),
            pixel_bounds,
        )
    }
}
