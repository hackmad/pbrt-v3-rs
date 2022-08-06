//! Direct Lighting Integrator

#![allow(dead_code)]

use atom::Atom;
use core::app::OPTIONS;
use core::camera::*;
use core::geometry::*;
use core::image_io::write_image;
use core::integrator::*;
use core::interaction::*;
use core::light::*;
use core::low_discrepency::radical_inverse;
use core::material::*;
use core::parallel::AtomicFloat;
use core::paramset::*;
use core::pbrt::*;
use core::reflection::*;
use core::sampler::*;
use core::scene::*;
use core::spectrum::*;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use samplers::HaltonSampler;
use std::sync::atomic::{AtomicI32, Ordering};
use std::sync::Arc;

/// Implements the stochastic progressive photon mapping integrator.
pub struct SPPMIntegrator {
    /// The camera.
    camera: ArcCamera,

    /// Initial search radius for photons.
    initial_search_radius: Float,

    /// Number of iterations of integrations to perform.
    n_iterations: usize,

    /// Maximum recursion depth.
    max_depth: usize,

    /// The total number of photons to be traced for the current iteration.
    photons_per_iteration: usize,

    /// How often SPPM image is stored in film.
    write_frequency: usize,
}

impl SPPMIntegrator {
    /// Create a new `SPPMIntegrator`.
    ///
    /// * `camera`                - The camera.
    /// * `initial_search_radius` - Initial search radius for photons.
    /// * `n_iterations`          - Number of iterations of integrations to perform.
    /// * `max_depth`             - Maximum recursion depth.
    /// * `photons_per_iteration` - The total number of photons to be traced for
    ///                             the current iteration.
    /// * `write_frequency`       - How often SPPM image is stored in film.
    pub fn new(
        camera: ArcCamera,
        initial_search_radius: Float,
        n_iterations: usize,
        max_depth: usize,
        photons_per_iteration: usize,
        write_frequency: usize,
    ) -> Self {
        Self {
            camera,
            initial_search_radius,
            n_iterations,
            max_depth,
            photons_per_iteration,
            write_frequency,
        }
    }
}

impl Integrator for SPPMIntegrator {
    /// Render the scene.
    ///
    /// * `scene` - The scene.
    fn render(&mut self, scene: &Scene) {
        // Initialize _pixelBounds_ and _pixels_ array for SPPM.
        let camera_data = self.camera.get_data();
        let pixel_bounds = camera_data.film.cropped_pixel_bounds;
        let n_pixels = pixel_bounds.area() as usize;

        let inv_sqrt_spp = 1.0 / (self.n_iterations as Float).sqrt();

        // Compute _lightDistr_ for sampling lights proportional to power.
        let light_distr = compute_light_power_distribution(scene);

        // Perform _nIterations_ of SPPM integration
        let sampler = HaltonSampler::new(self.n_iterations, pixel_bounds, false);

        // Compute number of tiles to use for SPPM camera pass.
        let pixel_extent = pixel_bounds.diagonal();
        let tile_size = OPTIONS.tile_size;
        let n_tiles = Point2::new(
            ((pixel_extent.x + tile_size as Int - 1) / tile_size as Int) as usize,
            ((pixel_extent.y + tile_size as Int - 1) / tile_size as Int) as usize,
        );
        let tile_count = n_tiles.x * n_tiles.y;

        let progress = if OPTIONS.quiet {
            ProgressBar::hidden()
        } else {
            let progress_style = ProgressStyle::default_bar()
                .template(
                    "{msg:25.cyan.bold} [{bar:40.green/white}] {pos:>5}/{len:5} ({elapsed}|{eta})",
                )
                .progress_chars("█▓▒░  ");
            let pb = ProgressBar::new(2 * self.n_iterations as u64);
            pb.set_message("Rendering");
            pb.set_style(progress_style);
            pb
        };

        let mut tiles: Vec<SPPMTile> = (0..tile_count)
            .map(|tile_idx| {
                // Render section of image corresponding to `tile`.
                let tile_x = tile_idx % n_tiles.x as usize;
                let tile_y = tile_idx / n_tiles.x as usize;

                // Compute _tileBounds_ for SPPM tile
                let x0 = pixel_bounds.p_min.x + tile_x as Int * tile_size as Int;
                let x1 = min(x0 + tile_size as Int, pixel_bounds.p_max.x);
                let y0 = pixel_bounds.p_min.y + tile_y as Int * tile_size as Int;
                let y1 = min(y0 + tile_size as Int, pixel_bounds.p_max.y);
                let tile_bounds = Bounds2i::new(Point2i::new(x0, y0), Point2i::new(x1, y1));

                let p_pixel_o_min = Point2i::from(tile_bounds.p_min - pixel_bounds.p_min);
                let p_pixel_o_max = Point2i::from(tile_bounds.p_max - pixel_bounds.p_min);
                let pixel_offset_min = (p_pixel_o_min.x
                    + p_pixel_o_min.y * (pixel_bounds.p_max.x - pixel_bounds.p_min.x))
                    as usize;
                let pixel_offset_max = (p_pixel_o_max.x
                    + p_pixel_o_min.y * (pixel_bounds.p_max.x - pixel_bounds.p_min.x))
                    as usize;

                SPPMTile::new(
                    (pixel_offset_min..pixel_offset_max) /* TODO double check if we are off by one / overlaps */
                        .map(|_| SPPMPixel::new(self.initial_search_radius))
                        .collect(),
                    tile_x,
                    tile_y,
                    pixel_offset_min,
                    Sampler::clone(&sampler, tile_idx as u64),
                    tile_bounds,
                )
            })
            .collect();

        for iter in 0..self.n_iterations {
            // Generate SPPM visible points
            tiles.par_iter_mut().for_each(|tile| {
                // Render section of image corresponding to `tile`.
                progress.set_message(format!("Rendering tile ({},{})", tile.x, tile.y));

                // Follow camera paths for _tile_ in image for SPPM.
                for p_pixel in tile.bounds {
                    let camera_sample = {
                        // Prepare `tile_sampler` for `p_pixel`.
                        let tile_sampler = Arc::get_mut(&mut tile.sampler).unwrap();
                        tile_sampler.start_pixel(&p_pixel);
                        tile_sampler.set_sample_number(iter);

                        // Generate camera ray for pixel for SPPM.
                        tile_sampler.get_camera_sample(&p_pixel)
                    };

                    let (mut ray, beta) = self.camera.generate_ray_differential(&camera_sample);
                    let mut beta = Spectrum::new(beta);
                    if beta.is_black() {
                        continue;
                    }
                    ray.scale_differentials(inv_sqrt_spp);

                    // Follow camera ray path until a visible point is created.

                    // Get `SPPMPixel` for `p_pixel`.
                    let p_pixel_o = Point2i::from(p_pixel - pixel_bounds.p_min);
                    let pixel_offset = (p_pixel_o.x
                        + p_pixel_o.y * (pixel_bounds.p_max.x - pixel_bounds.p_min.x))
                        as usize
                        - tile.offset_min; // Need to subtract tile.offset_min since slice starts at 0.

                    let mut specular_bounce = false;
                    let mut depth = 0;
                    while depth < self.max_depth {
                        let isect = scene.intersect(&mut ray);
                        if isect.is_none() {
                            // Accumulate light contributions for ray with no
                            // intersection.
                            for light in scene.lights.iter() {
                                tile.pixels[pixel_offset].ld += beta * light.le(&ray);
                            }
                            break;
                        }
                        let mut isect = isect.unwrap();

                        // Process SPPM camera ray intersection.

                        // Compute BSDF at SPPM camera ray intersection.
                        let mut bsdf: Option<BSDF> = None;
                        let mut bssrdf: Option<BSDF> = None;
                        isect.compute_scattering_functions(
                            &mut ray,
                            true,
                            TransportMode::Radiance,
                            &mut bsdf,
                            &mut bssrdf,
                        );
                        if bsdf.is_none() {
                            ray = isect.spawn_ray(&ray.d);
                            depth -= 1;
                            continue;
                        }
                        let bsdf = bsdf.unwrap();

                        // Accumulate direct illumination at SPPM camera ray
                        // intersection.
                        let wo = -ray.d;
                        if depth == 0 || specular_bounce {
                            tile.pixels[pixel_offset].ld += beta * isect.le(&wo);
                        }

                        let it = Interaction::Surface { si: isect };

                        tile.pixels[pixel_offset].ld += beta
                            * uniform_sample_one_light(
                                &it,
                                Some(&bsdf),
                                scene,
                                &mut tile.sampler,
                                false,
                                None,
                            );

                        // Need to move the `isect` back out of `it`.
                        let isect = match it {
                            Interaction::Surface { si } => si,
                            _ => unreachable!(), // isect was SurfaceInteraction. So this should not be possible.
                        };

                        // Possibly create visible point and end camera path.
                        let is_diffuse = bsdf.num_components(
                            BxDFType::BSDF_DIFFUSE
                                | BxDFType::BSDF_REFLECTION
                                | BxDFType::BSDF_TRANSMISSION,
                        ) > 0;
                        let is_glossy = bsdf.num_components(
                            BxDFType::BSDF_GLOSSY
                                | BxDFType::BSDF_REFLECTION
                                | BxDFType::BSDF_TRANSMISSION,
                        ) > 0;
                        if is_diffuse || (is_glossy && depth == self.max_depth - 1) {
                            tile.pixels[pixel_offset].vp =
                                VisiblePoint::new(isect.hit.p, wo, Some(bsdf), beta);
                            break;
                        }

                        // Spawn ray from SPPM camera path vertex.
                        if depth < self.max_depth - 1 {
                            let tile_sampler = Arc::get_mut(&mut tile.sampler).unwrap();
                            let BxDFSample {
                                f,
                                pdf,
                                wi,
                                bxdf_type,
                            } = bsdf.sample_f(&wo, &tile_sampler.get_2d(), BxDFType::all());
                            if pdf == 0.0 || f.is_black() {
                                break;
                            }

                            specular_bounce =
                                (bxdf_type & BxDFType::BSDF_SPECULAR) > BxDFType::BSDF_NONE;

                            beta *= f * wi.abs_dot(&isect.shading.n) / pdf;
                            if beta.y() < 0.25 {
                                let continue_prob = min(1.0, beta.y());
                                if tile_sampler.get_1d() > continue_prob {
                                    break;
                                }
                                beta /= continue_prob;
                            }
                            ray = isect.spawn_ray(&wi);
                        }

                        depth += 1;
                    }
                }
            });
            progress.inc(1);

            let mut grid_res: [Int; 3] = [0; 3];
            let mut grid_bounds = Bounds3f::default();

            // Allocate grid for SPPM visible points.
            let hash_size = n_pixels;
            let grid: Vec<Atom<Arc<SPPMPixelListNode>>> =
                (0..hash_size).map(|_| Atom::empty()).collect();

            // Compute grid bounds for SPPM visible points.
            let mut max_radius = 0.0;
            for tile in tiles.iter() {
                for pixel in tile.pixels.iter() {
                    if pixel.vp.beta.is_black() {
                        continue;
                    }
                    let vp_bound = Bounds3f::from(pixel.vp.p).expand(pixel.radius);
                    grid_bounds = grid_bounds.union(&vp_bound);
                    max_radius = max(max_radius, pixel.radius);
                }
            }

            // Compute resolution of SPPM grid in each dimension.
            let diag = grid_bounds.diagonal();
            let max_diag = diag.max_component();
            let base_grid_res = (max_diag / max_radius) as Int;
            assert!(base_grid_res > 0);
            for i in 0..3 {
                grid_res[i] = max((base_grid_res as Float * diag[i] / max_diag) as Int, 1);
            }

            // Add visible points to SPPM grid.
            tiles.par_iter().chunks(4096 / tile_size).for_each(|tiles| {
                for tile in tiles {
                    for pixel in tile.pixels.iter() {
                        if !pixel.vp.beta.is_black() {
                            // Add pixel's visible point to applicable grid cells.
                            let radius = pixel.radius;
                            let (_, p_min) = to_grid(
                                &(pixel.vp.p - Vector3f::new(radius, radius, radius)),
                                &grid_bounds,
                                &grid_res,
                            );
                            let (_, p_max) = to_grid(
                                &(pixel.vp.p + Vector3f::new(radius, radius, radius)),
                                &grid_bounds,
                                &grid_res,
                            );
                            for z in p_min.z..=p_max.z {
                                for y in p_min.y..=p_max.y {
                                    for x in p_min.x..=p_max.x {
                                        // Add visible point to grid cell `(x, y, z)`.
                                        let h = hash(&Point3i::new(x, y, z), hash_size);
                                        let mut node = Arc::new(SPPMPixelListNode::new(&pixel));

                                        // Atomically add `node` to the start of `grid[h]`'s linked list.
                                        let list =
                                            grid[h].swap(Arc::clone(&node), Ordering::AcqRel);
                                        if let Some(list) = list {
                                            Arc::get_mut(&mut node).unwrap().next = Some(list);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            });

            // Trace photons and accumulate contributions.
            (0..self.photons_per_iteration)
                .into_par_iter()
                .chunks(8192)
                .for_each(|chunk| {
                    for photon_index in chunk {
                        // Follow photon path for `photon_index`.
                        let halton_index =
                            (iter * self.photons_per_iteration + photon_index) as u64;
                        let mut halton_dim = 0;

                        // Choose light to shoot photon from.
                        let light_sample = radical_inverse(halton_dim, halton_index);
                        halton_dim += 1;

                        let mut beta = Spectrum::ZERO;
                        let mut photon_ray = Ray::default();
                        if let Some(distr) = light_distr.as_ref() {
                            let (light_num, light_pdf, _u_remapped) =
                                distr.sample_discrete(light_sample);
                            let light = &scene.lights[light_num];

                            // Compute sample values for photon ray leaving light source.
                            let u_light_0 = Point2f::new(
                                radical_inverse(halton_dim, halton_index),
                                radical_inverse(halton_dim + 1, halton_index),
                            );
                            let u_light_1 = Point2f::new(
                                radical_inverse(halton_dim + 2, halton_index),
                                radical_inverse(halton_dim + 3, halton_index),
                            );
                            let u_light_time = lerp(
                                radical_inverse(halton_dim + 4, halton_index),
                                camera_data.shutter_open,
                                camera_data.shutter_close,
                            );
                            halton_dim += 5;

                            // Generate _photonRay_ from light source and initialize _beta_
                            let Le {
                                ray,
                                n_light,
                                pdf_pos,
                                pdf_dir,
                                value: le,
                            } = light.sample_le(&u_light_0, &u_light_1, u_light_time);
                            if pdf_pos == 0.0 || pdf_dir == 0.0 || le.is_black() {
                                return;
                            }

                            beta = (n_light.abs_dot(&ray.d) * le) / (light_pdf * pdf_pos * pdf_dir);
                            photon_ray = ray;
                        }
                        if beta.is_black() {
                            return;
                        }

                        // Follow photon path through scene and record intersections.
                        let mut depth = 0;
                        while depth < self.max_depth {
                            let isect = scene.intersect(&mut photon_ray);
                            if isect.is_none() {
                                break;
                            }
                            let mut isect = isect.unwrap();

                            if depth > 0 {
                                // Add photon contribution to nearby visible points.
                                let (result, photon_grid_index) =
                                    to_grid(&isect.hit.p, &grid_bounds, &grid_res);
                                if result {
                                    let h = hash(&photon_grid_index, hash_size);

                                    // Add photon contribution to visible points in `grid[h]`.
                                    let mut current_node = grid[h].take(Ordering::Relaxed);
                                    while let Some(node) = current_node.as_ref() {
                                        let pixel = node.pixel;
                                        let radius = pixel.radius;
                                        if pixel.vp.p.distance_squared(isect.hit.p)
                                            > radius * radius
                                        {
                                            continue;
                                        }
                                        // Update `pixel` `Phi` and `M` for nearby photon.
                                        let wi = -photon_ray.d;
                                        let phi = beta
                                            * pixel
                                                .vp
                                                .bsdf
                                                .as_ref()
                                                .map_or(Spectrum::ZERO, |bsdf| {
                                                    bsdf.f(&pixel.vp.wo, &wi, BxDFType::all())
                                                });
                                        for i in 0..SPECTRUM_SAMPLES {
                                            pixel.phi[i].add(phi[i]);
                                        }
                                        pixel.m.fetch_add(1, Ordering::Relaxed);

                                        current_node = node.next.as_ref().map(Arc::clone);
                                    }
                                }
                            }
                            // Sample new photon ray direction.

                            // Compute BSDF at photon intersection point.
                            let mut bsdf: Option<BSDF> = None;
                            let mut bssrdf: Option<BSDF> = None;
                            isect.compute_scattering_functions(
                                &mut photon_ray,
                                true,
                                TransportMode::Importance,
                                &mut bsdf,
                                &mut bssrdf,
                            );
                            if bsdf.is_none() {
                                depth -= 1;
                                photon_ray = isect.spawn_ray(&photon_ray.d);
                                continue;
                            }
                            let photon_bsdf = bsdf.unwrap();

                            // Sample BSDF `fr` and direction `wi` for reflected photon.
                            let wo = -photon_ray.d;

                            // Generate `bsdf_sample` for outgoing photon sample.
                            let bsdf_sample = Point2f::new(
                                radical_inverse(halton_dim, halton_index),
                                radical_inverse(halton_dim + 1, halton_index),
                            );
                            halton_dim += 2;
                            let BxDFSample {
                                f: fr,
                                pdf,
                                wi,
                                bxdf_type: _flags,
                            } = photon_bsdf.sample_f(&wo, &bsdf_sample, BxDFType::all());
                            if fr.is_black() || pdf == 0.0 {
                                break;
                            }
                            let bnew = beta * fr * wi.abs_dot(&isect.shading.n) / pdf;

                            // Possibly terminate photon path with Russian roulette.
                            let q = max(0.0, 1.0 - bnew.y() / beta.y());
                            let rad_inv = radical_inverse(halton_dim, halton_index);
                            halton_dim += 1;
                            if rad_inv < q {
                                break;
                            }
                            beta = bnew / (1.0 - q);
                            photon_ray = isect.spawn_ray(&wi);

                            depth += 1;
                        }
                    }
                });
            progress.inc(1);

            // Update pixel values from this pass' photons.
            tiles
                .par_iter_mut()
                .chunks(4096 / tile_size)
                .for_each(|tiles| {
                    for tile in tiles {
                        for p in tile.pixels.iter() {
                            let m = p.m.load(Ordering::Relaxed);
                            if m > 0 {
                                // Update pixel photon count, search radius, and `tau` from
                                // photons.
                                const GAMMA: Float = 2.0 / 3.0;
                                let n_new = p.n + GAMMA * m as Float;
                                let r_new = p.radius * (n_new / (p.n + m as Float)).sqrt();

                                let mut phi = Spectrum::ZERO;
                                for j in 0..SPECTRUM_SAMPLES {
                                    phi[j] = p.phi[j].load(Ordering::Relaxed);
                                }

                                p.tau = (p.tau + p.vp.beta * phi) * (r_new * r_new)
                                    / (p.radius * p.radius);
                                p.n = n_new;
                                p.radius = r_new;
                                p.m.store(0, Ordering::Relaxed);

                                for j in 0..SPECTRUM_SAMPLES {
                                    p.phi[j].store(0.0, Ordering::Relaxed);
                                }
                            }
                            // Reset _VisiblePoint_ in pixel
                            p.vp.beta = Spectrum::ZERO;
                            p.vp.bsdf = None;
                        }
                    }
                });

            // Periodically store SPPM image in film and write image.
            if iter + 1 == self.n_iterations || (iter + 1) % self.write_frequency == 0 {
                let x0 = pixel_bounds.p_min.x;
                let x1 = pixel_bounds.p_max.x;
                let Np = (iter + 1) as u64 * self.photons_per_iteration as u64;
                let image = vec![Spectrum::ZERO; n_pixels];
                let mut offset = 0;
                for y in pixel_bounds.p_min.y..pixel_bounds.p_max.y {
                    for x in x0..x1 {
                        // Compute radiance `L` for SPPM pixel `pixel`.
                        let pixel = &pixels[(y - pixel_bounds.p_min.y) * (x1 - x0) + (x - x0)];
                        let mut L = pixel.ld / (iter + 1);
                        L += pixel.tau / (Np * PI * pixel.radius * pixel.radius);
                        image[offset] = L;
                        offset += 1;
                    }
                }
                camera_data.film.set_image(&image);
                camera_data.film.write_image(1.0);
                // Write SPPM radius image, if requested.

                if OPTIONS.sppm_radius {
                    let mut rimg = vec![0.0; 3 * n_pixels];
                    let minrad = 1e30;
                    let maxrad = 0.0;
                    for y in pixel_bounds.p_min.y..pixel_bounds.p_max.y {
                        for x in x0..x1 {
                            let p = &pixels[(y - pixel_bounds.p_min.y) * (x1 - x0) + (x - x0)];
                            minrad = min(minrad, p.radius);
                            maxrad = max(maxrad, p.radius);
                        }
                    }
                    debug!(
                        "Iterations: {} ({:2} s) radius range: {} - {}",
                        iter + 1,
                        progress.elapsed().as_secs_f32(),
                        minrad,
                        maxrad
                    );
                    let mut offset = 0;
                    for y in pixel_bounds.p_min.y..pixel_bounds.p_max.y {
                        for x in x0..x1 {
                            let p = &pixels[(y - pixel_bounds.p_min.y) * (x1 - x0) + (x - x0)];
                            let v = 1.0 - (p.radius - minrad) / (maxrad - minrad);
                            rimg[offset] = v;
                            rimg[offset + 1] = v;
                            rimg[offset + 2] = v;
                            offset += 3;
                        }
                    }
                    if let Err(err) = write_image("sppm_radius.png", &rimg, &pixel_bounds) {
                        error!("Error writing output image sppm_radius.png. {:}.", err);
                    }
                }
            }
        }

        progress.finish_with_message("Done");
    }

    fn li(
        &self,
        _ray: &mut Ray,
        _scene: &Scene,
        _sampler: &mut ArcSampler,
        _depth: usize,
    ) -> Spectrum {
        Spectrum::ZERO
    }
}

impl From<(&ParamSet, ArcCamera)> for SPPMIntegrator {
    /// Create a `SPPMIntegrator` from given parameter set and camera.
    ///
    /// * `p` - A tuple containing parameter set and camera.
    fn from(p: (&ParamSet, ArcCamera)) -> Self {
        let (params, camera) = p;

        let n_iterations = params.find_one_int("numiterations", 64);
        let mut n_iterations = params.find_one_int("iterations", n_iterations) as usize;
        if OPTIONS.quick_render {
            n_iterations = max(1, n_iterations / 16);
        }

        let max_depth = params.find_one_int("maxdepth", 5) as usize;
        let photons_per_iteration = params.find_one_int("photonsperiteration", -1) as usize;
        let write_frequency = params.find_one_int("imagewritefrequency", 1 << 31) as usize;
        let initial_search_radius = params.find_one_float("radius", 1.0);

        Self::new(
            camera,
            initial_search_radius,
            n_iterations,
            max_depth,
            photons_per_iteration,
            write_frequency,
        )
    }
}

struct SPPMTile {
    pixels: Vec<SPPMPixel>,
    x: usize,
    y: usize,
    offset_min: usize,
    sampler: ArcSampler,
    bounds: Bounds2i,
}

impl SPPMTile {
    fn new(
        pixels: Vec<SPPMPixel>,
        tile_x: usize,
        tile_y: usize,
        offset_min: usize,
        sampler: ArcSampler,
        bounds: Bounds2i,
    ) -> Self {
        Self {
            pixels,
            x: tile_x,
            y: tile_y,
            offset_min,
            sampler,
            bounds,
        }
    }
}

#[derive(Default)]
struct SPPMPixel {
    radius: Float,
    ld: Spectrum,
    vp: VisiblePoint,
    phi: [AtomicFloat; SPECTRUM_SAMPLES],
    m: AtomicI32,
    n: Float,
    tau: Spectrum,
}

impl SPPMPixel {
    /// Create a new `SPPMPixel` with initial search radius and default field
    /// values.
    fn new(radius: Float) -> Self {
        Self {
            radius,
            ..Default::default()
        }
    }
}

#[derive(Default, Clone)]
struct VisiblePoint {
    p: Point3f,
    wo: Vector3f,
    bsdf: Option<BSDF>,
    beta: Spectrum,
}

impl VisiblePoint {
    fn new(p: Point3f, wo: Vector3f, bsdf: Option<BSDF>, beta: Spectrum) -> Self {
        Self { p, wo, bsdf, beta }
    }
}

struct SPPMPixelListNode<'a> {
    pixel: &'a SPPMPixel,
    next: Option<Arc<SPPMPixelListNode<'a>>>,
}

impl<'a> SPPMPixelListNode<'a> {
    fn new(pixel: &'a SPPMPixel) -> Self {
        Self { pixel, next: None }
    }
}

fn to_grid(p: &Point3f, bounds: &Bounds3f, grid_res: &[Int; 3]) -> (bool, Point3i) {
    let mut pi = Point3i::default();
    let mut in_bounds = true;
    let pg = bounds.offset(&p);
    for i in 0..3 {
        pi[i] = (grid_res[i] as Float * pg[i]) as Int;
        in_bounds &= pi[i] >= 0 && pi[i] < grid_res[i];
        pi[i] = clamp(pi[i], 0, grid_res[i] - 1);
    }
    (in_bounds, pi)
}

#[inline(always)]
fn hash(p: &Point3i, hash_size: usize) -> usize {
    ((p.x * 73856093) ^ (p.y * 19349663) ^ (p.z * 83492791)) as usize % hash_size
}
