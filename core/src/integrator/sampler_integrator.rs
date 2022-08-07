//! Sampler Integrator

use super::*;
use crate::app::OPTIONS;
use crate::camera::*;
use crate::film::FilmTile;
use crate::geometry::*;
use crate::interaction::*;
use crate::pbrt::*;
use crate::reflection::*;
use crate::sampler::*;
use crate::scene::*;
use crate::spectrum::*;
use std::sync::Arc;

/// Common data for sampler integrators.
pub struct SamplerIntegratorData {
    /// Sampler responsible for choosing points on the image plane from which
    /// to trace rays and for supplying sample positions used by integrators.
    pub sampler: ArcSampler,

    /// The camera.
    pub camera: ArcCamera,

    /// Pixel bounds for the image.
    pub pixel_bounds: Bounds2i,

    /// Maximum recursion depth.
    pub max_depth: usize,
}

impl SamplerIntegratorData {
    /// Create a new `SamplerIntegratorData`.
    ///
    /// * `max_depth`    - Maximum recursion depth.
    /// * `camera`       - The camera.
    /// * `sampler`      - Sampler responsible for choosing point on image plane
    ///                    from which to trace rays.
    /// * `pixel_bounds` - Pixel bounds for the image.
    pub fn new(
        max_depth: usize,
        camera: ArcCamera,
        sampler: ArcSampler,
        pixel_bounds: Bounds2i,
    ) -> Self {
        Self {
            camera: Arc::clone(&camera),
            max_depth,
            sampler,
            pixel_bounds,
        }
    }
}

/// Implements the basis of a rendering process driven by a stream of samples
/// from a `Sampler`. Each sample identifies a point on the image plane at
/// which we compute the light arriving from the scene.
pub trait SamplerIntegrator: Integrator + Send + Sync {
    /// Returns the common data.
    fn get_data(&self) -> &SamplerIntegratorData;

    /// Trace rays for specular reflection.
    ///
    /// * `ray`     - The ray.
    /// * `isect`   - The surface interaction.
    /// * `bsdf`    - The BSDF at the surface interaction.
    /// * `scene`   - The scene.
    /// * `sampler` - The sampler.
    /// * `depth`   - The recursive depth.
    fn specular_reflect(
        &self,
        ray: &mut Ray,
        isect: &SurfaceInteraction,
        bsdf: &Option<BSDF>,
        scene: &Scene,
        sampler: &mut ArcSampler,
        depth: usize,
    ) -> Spectrum {
        if let Some(bsdf) = bsdf.as_ref() {
            // Compute specular reflection direction `wi` and BSDF value.
            let wo = isect.hit.wo;

            let sample = Arc::get_mut(sampler).unwrap().get_2d();
            let bxdf_type = BxDFType::BSDF_REFLECTION | BxDFType::BSDF_SPECULAR;
            let BxDFSample {
                f,
                pdf,
                wi,
                bxdf_type: _,
            } = bsdf.sample_f(&wo, &sample, bxdf_type);

            // Return contribution of specular reflection
            let ns = isect.shading.n;
            if pdf > 0.0 && !f.is_black() && wi.abs_dot(&ns) != 0.0 {
                // Compute ray differential `rd` for specular reflection.
                let mut rd = isect.hit.spawn_ray(&wi);
                if let Some(differentials) = ray.differentials {
                    let rx_origin = isect.hit.p + isect.der.dpdx;
                    let ry_origin = isect.hit.p + isect.der.dpdy;

                    // Compute differential reflected directions.
                    let dndx =
                        isect.shading.dndu * isect.der.dudx + isect.shading.dndv * isect.der.dvdx;
                    let dndy =
                        isect.shading.dndu * isect.der.dudy + isect.shading.dndv * isect.der.dvdy;
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

                return f * self.li(&mut rd, scene, sampler, depth + 1) * wi.abs_dot(&ns) / pdf;
            }
        }

        Spectrum::new(0.0)
    }

    /// Trace rays for specular refraction.
    ///
    /// * `ray`     - The ray.
    /// * `isect`   - The surface interaction.
    /// * `bsdf`    - The BSDF at the surface interaction.
    /// * `scene`   - The scene.
    /// * `sampler` - The sampler.
    /// * `depth`   - The recursive depth.
    fn specular_transmit(
        &self,
        ray: &mut Ray,
        isect: &SurfaceInteraction,
        bsdf: &Option<BSDF>,
        scene: &Scene,
        sampler: &mut ArcSampler,
        depth: usize,
    ) -> Spectrum {
        if let Some(bsdf) = bsdf.as_ref() {
            let wo = isect.hit.wo;
            let p = isect.hit.p;

            let sample = Arc::get_mut(sampler).unwrap().get_2d();
            let bxdf_type = BxDFType::BSDF_TRANSMISSION | BxDFType::BSDF_SPECULAR;
            let BxDFSample {
                f,
                pdf,
                wi,
                bxdf_type: _,
            } = bsdf.sample_f(&wo, &sample, bxdf_type);

            let mut ns = isect.shading.n;
            if pdf > 0.0 && !f.is_black() && wi.abs_dot(&ns) != 0.0 {
                // Compute ray differential _rd_ for specular transmission
                let mut rd = isect.hit.spawn_ray(&wi);
                if let Some(differentials) = ray.differentials {
                    let rx_origin = p + isect.der.dpdx;
                    let ry_origin = p + isect.der.dpdy;

                    let mut dndx =
                        isect.shading.dndu * isect.der.dudx + isect.shading.dndv * isect.der.dvdx;
                    let mut dndy =
                        isect.shading.dndu * isect.der.dudy + isect.shading.dndv * isect.der.dvdy;

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

                return f * self.li(&mut rd, scene, sampler, depth + 1) * wi.abs_dot(&ns) / pdf;
            }
        }

        Spectrum::new(0.0)
    }

    /// Render the scene.
    ///
    /// * `scene` - The scene.
    fn render(&self, scene: &Scene) {
        // Get sampler and camera data.
        let data = self.get_data();
        let camera = Arc::clone(&data.camera);
        let camera_data = camera.get_data();

        // Compute number of tiles, `n_tiles`, to use for parallel rendering.
        let sample_bounds = camera_data.film.get_sample_bounds();
        let sample_extent = sample_bounds.diagonal();
        let tile_size = OPTIONS.tile_size as i32;
        let n_tiles = Point2::new(
            ((sample_extent.x + tile_size - 1) / tile_size) as usize,
            ((sample_extent.y + tile_size - 1) / tile_size) as usize,
        );
        let tile_count = n_tiles.x * n_tiles.y;

        // Parallelize.
        info!("Rendering {}x{} tiles", n_tiles.x, n_tiles.y);

        let progress = create_progress_reporter(tile_count as u64 + 1_u64); // Render + image write
        progress.set_message("Rendering scene");

        crossbeam::scope(|scope| {
            let (tx, rx) = crossbeam_channel::bounded(OPTIONS.threads());

            // Spawn worker threads.
            for _ in 0..OPTIONS.threads() {
                let rxc = rx.clone();
                let progress = &progress;
                scope.spawn(move |_| {
                    for tile_idx in rxc.iter() {
                        // Render section of image corresponding to `tile`.
                        let film_tile = self.render_tile(tile_idx, n_tiles, scene, sample_bounds);

                        // Merge image tile into `Film`.
                        camera_data.film.merge_film_tile(&film_tile);
                        progress.inc(1);
                    }
                });
            }
            drop(rx); // Drop extra rx since we've cloned one for each woker.

            // Send work.
            for tile_idx in 0..tile_count {
                tx.send(tile_idx).unwrap();
            }
        })
        .unwrap();

        // Save final image after rendering.
        progress.set_message("Writing image");
        progress.inc(1);
        camera_data.film.write_image(1.0);

        progress.finish_with_message("Render complete");
    }

    /// Render an image tile.
    ///
    /// * `tile_idx`      - Unique tile index.
    /// * `n_tiles`       - Number of tiles in (x, y) direction.
    /// * `scene`         - Scene.
    /// * `sample_bounds` - Sample bounds.
    fn render_tile(
        &self,
        tile_idx: usize,
        n_tiles: Point2<usize>,
        scene: &Scene,
        sample_bounds: Bounds2i,
    ) -> FilmTile {
        // Get the x and y tile indices.
        let tile_x = tile_idx % n_tiles.x;
        let tile_y = tile_idx / n_tiles.x;

        // Get camera data.
        let data = self.get_data();
        let camera = Arc::clone(&data.camera);
        let camera_data = camera.get_data();

        // Get sampler instance for tile.
        let mut tile_sampler = Sampler::clone(&*data.sampler, tile_idx as u64);

        let samples_per_pixel = {
            let tile_sampler_data = Arc::get_mut(&mut tile_sampler).unwrap().get_data();
            tile_sampler_data.samples_per_pixel
        };

        // Compute sample bounds for tile.
        let tile_size = OPTIONS.tile_size as i32;
        let x0 = sample_bounds.p_min.x + tile_x as i32 * tile_size;
        let x1 = min(x0 + tile_size, sample_bounds.p_max.x);
        let y0 = sample_bounds.p_min.y + tile_y as i32 * tile_size;
        let y1 = min(y0 + tile_size, sample_bounds.p_max.y);
        let tile_bounds = Bounds2i::new(Point2i::new(x0, y0), Point2i::new(x1, y1));

        info!("Starting image tile ({tile_x}, {tile_y}) -> {tile_bounds}");

        let mut film_tile = camera_data.film.get_film_tile(tile_bounds);

        // Loop over pixels in tile to render them.
        for pixel in tile_bounds {
            Arc::get_mut(&mut tile_sampler).unwrap().start_pixel(&pixel);

            // Do this check after the StartPixel() call; this keeps the
            // usage of RNG values from (most) Samplers that use RNGs
            // consistent, which improves reproducability / debugging.
            if !data.pixel_bounds.contains_exclusive(&pixel) {
                continue;
            }

            loop {
                // Initialize `CameraSample` for current sample.
                let camera_sample = Arc::get_mut(&mut tile_sampler)
                    .unwrap()
                    .get_camera_sample(&pixel);

                // Generate camera ray for current sample.
                let (mut ray, ray_weight) = { camera.generate_ray_differential(&camera_sample) };
                ray.scale_differentials(1.0 / (samples_per_pixel as Float).sqrt());

                // Evaluate radiance along camera ray.
                let mut l = Spectrum::new(0.0);
                if ray_weight > 0.0 {
                    l = self.li(&mut ray, scene, &mut tile_sampler, 0);
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

                debug!("Pixel: {pixel}, Camera sample: {camera_sample} -> ray: {ray}, ray weight {ray_weight} -> L = {l}");

                // Add camera ray's contribution to image.
                film_tile.add_sample(camera_sample.p_film, l, ray_weight);

                if !Arc::get_mut(&mut tile_sampler).unwrap().start_next_sample() {
                    break;
                }
            }
        }

        info!("Finished image tile ({tile_x}, {tile_y}) -> {tile_bounds}");

        film_tile
    }
}
