//! Sampler Integrator

use crate::core::camera::*;
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::reflection::*;
use crate::core::sampler::*;
use crate::core::scene::*;
use crate::core::spectrum::*;
use rayon::prelude::*;
use std::sync::Arc;
use itertools::iproduct;

/// Implements the basis of a rendering process driven by a stream of samples
/// from a `Sampler`. Each sample identifies a point on the image plane at
/// which we compute the light arriving from the scene.
pub struct SamplerIntegrator {
    /// Sampler responsible for choosing points on the image plane from which
    /// to trace rays and for supplying sample positions used by integrators.
    pub sampler: ArcSampler,

    /// The camera.
    pub camera: ArcCamera,

    /// Pixel bounds for the image.
    pub pixel_bounds: Bounds2i,
}

impl SamplerIntegrator {
    /// Create a new `SamplerIntegrator`.
    ///
    /// * `camera`       - The camera.
    /// * `sampler`      - Sampler responsible for choosing point on image plane
    ///                    from which to trace rays.
    /// * `pixel_bounds` - Pixel bounds for the image.
    pub fn new(camera: ArcCamera, sampler: ArcSampler, pixel_bounds: Bounds2i) -> Self {
        Self {
            camera,
            sampler,
            pixel_bounds,
        }
    }

    /// Trace rays for specular reflection.
    ///
    /// * `ray`   - The ray.
    /// * `isect` - The surface interaction.
    /// * `scene` - The scene.
    /// * `depth` - The recursive depth.
    /// * `li`    - A lambda that returns the incident radiance at the origin of
    ///             a given ray.
    pub fn specular_reflect(
        &mut self,
        ray: &mut Ray,
        isect: &SurfaceInteraction,
        scene: Arc<Scene>,
        depth: usize,
        li: &dyn Fn(&Ray, Arc<Scene>, &mut ArcSampler, usize) -> Spectrum,
    ) -> Spectrum {
        if let Some(bsdf) = isect.bsdf.clone() {
            // Compute specular reflection direction `wi` and BSDF value.
            let wo = isect.hit.wo;

            let sample = Arc::get_mut(&mut self.sampler).unwrap().get_2d();
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
                return f * li(&rd, scene.clone(), &mut self.sampler, depth + 1) * wi.abs_dot(&ns)
                    / pdf;
            }
        }

        Spectrum::new(0.0)
    }

    /// Trace rays for specular refraction.
    ///
    /// * `ray`   - The ray.
    /// * `isect` - The surface interaction.
    /// * `scene` - The scene.
    /// * `depth` - The recursive depth.
    /// * `li`    - A lambda that returns the incident radiance at the origin of
    ///             a given ray.
    pub fn specular_transmit(
        &mut self,
        ray: &mut Ray,
        isect: &SurfaceInteraction,
        scene: Arc<Scene>,
        depth: usize,
        li: &dyn Fn(&Ray, Arc<Scene>, &mut ArcSampler, usize) -> Spectrum,
    ) -> Spectrum {
        if let Some(bsdf) = isect.bsdf.clone() {
            let wo = isect.hit.wo;
            let p = isect.hit.p;

            let sample = Arc::get_mut(&mut self.sampler).unwrap().get_2d();
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
                return f * li(&rd, scene.clone(), &mut self.sampler, depth + 1) * wi.abs_dot(&ns)
                    / pdf;
            }
        }

        Spectrum::new(0.0)
    }
    
    /// Render the scene.
    /// 
    /// NOTE: The integrators that use this function should call their own
    /// preprocess(scene, sampler) implementation before calling this.
    ///
    /// * `scene` - The scene.
    /// * `li`    - A lambda that returns the incident radiance at the origin of
    ///             a given ray.
    pub fn render(
        &mut self,
        scene: Arc<Scene>, 
        li: &'static (dyn Fn(&Ray, Arc<Scene>, &mut ArcSampler, usize) -> Spectrum + Sync)) {
        // Compute number of tiles, `n_tiles`, to use for parallel rendering
        let film = self.camera.get_data().film.clone();
        let sample_bounds = film.get_sample_bounds();
        let sample_extent = sample_bounds.diagonal();
        let tile_size = 16;
        let n_tiles = Point2::new(
            ((sample_extent.x + tile_size - 1) / tile_size) as usize, 
            ((sample_extent.y + tile_size - 1) / tile_size) as usize
        );

        info!("Rendering tiles {}", n_tiles.x * n_tiles.y);
       
        // Parallelize.
        let tiles = iproduct!(0..n_tiles.x, 0..n_tiles.y).par_bridge();
        tiles.for_each(|t| {
            // Render section of image corresponding to `tile`.
            let tile = Point2::new(t.0, t.1);

            // Get sampler instance for tile.
            let seed = tile.y * n_tiles.x + tile.x;
            let mut tile_sampler = Sampler::clone(&*self.sampler, seed as u64);

            let samples_per_pixel = {
                let tile_sampler_data = Arc::get_mut(&mut tile_sampler).unwrap().get_data();
                tile_sampler_data.samples_per_pixel
            };

            // Compute sample bounds for tile
            let x0 = sample_bounds.p_min.x + tile.x as i32 * tile_size;
            let x1 = min(x0 + tile_size, sample_bounds.p_max.x);
            let y0 = sample_bounds.p_min.y + tile.y as i32 * tile_size;
            let y1 = min(y0 + tile_size, sample_bounds.p_max.y);
            let tile_bounds = Bounds2i::new(Point2i::new(x0, y0), Point2i::new(x1, y1));
            info!("Starting image tile {:}" , tile_bounds);

            // Get _FilmTile_ for tile
            let mut film_tile = film.get_film_tile(tile_bounds);

            // Loop over pixels in tile to render them
            for pixel in tile_bounds {
                Arc::get_mut(&mut tile_sampler).unwrap().start_pixel(&pixel);

                // Do this check after the StartPixel() call; this keeps
                // the usage of RNG values from (most) Samplers that use
                // RNGs consistent, which improves reproducability /
                // debugging.

                if !self.pixel_bounds.contains_exclusive(&pixel) {
                    continue;
                }

                loop {
                    // Initialize _CameraSample_ for current sample
                    let camera_sample = Arc::get_mut(&mut tile_sampler).unwrap().get_camera_sample(&pixel);

                    // Generate camera ray for current sample
                    let (mut ray, ray_weight) = self.camera.generate_ray_differential(&camera_sample);
                    ray.scale_differentials(1.0 / (samples_per_pixel as Float).sqrt());

                    // Evaluate radiance along camera ray
                    let mut l = Spectrum::new(0.0);
                    if ray_weight > 0.0 {
                        l = li(&ray, scene.clone(), &mut tile_sampler, 0);
                    }

                    // Issue warning if unexpected radiance value returned
                    let tile_sampler_data = Arc::get_mut(&mut tile_sampler).unwrap().get_data();
                    let current_sample_number = tile_sampler_data.current_sample_number();
                    if l.has_nans() {
                        error!("Not-a-number radiance value returned for pixel ({}, {}), sample {}. Setting to black.",
                            pixel.x, pixel.y, current_sample_number);
                        l = Spectrum::new(0.0);
                    } else if l.y() < -1e-5 {
                        error!("Negative luminance value, {}, returned for pixel ({}, {}), sample {}. Setting to black.",
                            l.y(), pixel.x, pixel.y, current_sample_number);
                        l = Spectrum::new(0.0);
                    } else if l.y().is_infinite() {
                        error!("Infinite luminance value returned for pixel ({}, {}), sample {}. Setting to black.",
                            pixel.x, pixel.y, current_sample_number);
                        l = Spectrum::new(0.0);
                    }

                    info!("Camera sample: {:} -> ray: {:} -> L = {:}", camera_sample, ray, l);

                    // Add camera ray's contribution to image
                    Arc::get_mut(&mut film_tile).unwrap().add_sample(camera_sample.p_film, l, ray_weight);

                    if !Arc::get_mut(&mut tile_sampler).unwrap().start_next_sample() {
                        break;
                    }
                } 
            }
            info!("Finished image tile {:}", tile_bounds);

            // Merge image tile into `Film`.
            let mut film_mut = film.clone();
            Arc::get_mut(&mut film_mut).unwrap().merge_film_tile(film_tile.clone());
        });
        
        info!("Rendering finished.");

        // Save final image after rendering.
        film.clone().write_image(1.0);
        info!("Output image written.");
    }
}
