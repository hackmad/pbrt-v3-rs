//! Direct Lighting Integrator

#![allow(dead_code)]

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
use core::sampling::Distribution1D;
use core::scene::*;
use core::spectrum::*;
use samplers::HaltonSampler;
use std::collections::VecDeque;
use std::sync::atomic::{AtomicI32, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;

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
        let photons_per_iter = if photons_per_iteration > 0 {
            photons_per_iteration
        } else {
            camera.get_data().film.cropped_pixel_bounds.area() as usize
        };
        Self {
            camera,
            initial_search_radius,
            n_iterations,
            max_depth,
            photons_per_iteration: photons_per_iter,
            write_frequency,
        }
    }
}

impl Integrator for SPPMIntegrator {
    /// Preprocess the scene.
    ///
    /// * `scene` - The scene
    fn preprocess(&mut self, _scene: &Scene) {}

    /// Render the scene.
    ///
    /// * `scene` - The scene.
    fn render(&self, scene: &Scene) {
        // Initialize _pixelBounds_ and _pixels_ array for SPPM.
        let camera_data = self.camera.get_data();
        let pixel_bounds = camera_data.film.cropped_pixel_bounds;
        let n_pixels = pixel_bounds.area() as usize;
        let mut pixels: Vec<SPPMPixel> = (0..n_pixels)
            .into_iter()
            .map(|_| SPPMPixel::new(self.initial_search_radius))
            .collect();
        let inv_sqrt_spp = 1.0 / (self.n_iterations as Float).sqrt();

        // Compute `light_distr` for sampling lights proportional to power.
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

        let progress = create_progress_bar(2 * self.n_iterations as u64);
        progress.set_message("Rendering scene");

        for iter in 0..self.n_iterations {
            // Generate SPPM visible points.
            generate_visible_points(
                &mut pixels,
                iter,
                &n_tiles,
                tile_size,
                tile_count,
                scene,
                &sampler,
                &pixel_bounds,
                inv_sqrt_spp,
                self.max_depth,
                Arc::clone(&self.camera),
            );
            progress.inc(1);

            // Allocate grid for SPPM visible points.
            let hash_size = n_pixels;
            let grid: Vec<Arc<Mutex<VecDeque<&SPPMPixel>>>> =
                (0..hash_size).map(|_| Arc::new(Mutex::new(VecDeque::new()))).collect();

            // Compute grid bounds for SPPM visible points.
            let (grid_bounds, max_radius) = compute_grid_bounds(&pixels);

            // Compute resolution of SPPM grid in each dimension.
            let grid_res = compute_grid_resolution(&grid_bounds, max_radius);

            // Add visible points to SPPM grid.
            add_visible_points_to_grid(&pixels, &grid_bounds, &grid_res, hash_size, &grid);

            // Trace photons and accumulate contributions.
            trace_photons(
                iter,
                self.photons_per_iteration,
                scene,
                camera_data,
                &light_distr,
                self.max_depth,
                &grid_bounds,
                &grid_res,
                hash_size,
                &grid,
            );
            progress.inc(1);

            // Update pixel values from this pass' photons.
            update_pixels(&mut pixels);

            // Periodically store SPPM image in film and write image.
            if iter + 1 == self.n_iterations || (iter + 1) % self.write_frequency == 0 {
                write_sppm_image(
                    &pixels,
                    iter,
                    &pixel_bounds,
                    self.photons_per_iteration,
                    &camera_data,
                    progress.elapsed(),
                );
            }
        }

        progress.finish_with_message("Render complete");
    }

    /// Returns the incident radiance at the origin of a given ray.
    ///
    /// NOTE: This is never called.
    ///
    /// * `ray`     - The ray.
    /// * `scene`   - The scene.
    /// * `sampler` - The sampler.
    /// * `depth`   - The recursion depth.
    fn li(&self, _ray: &mut Ray, _scene: &Scene, _sampler: &mut ArcSampler, _depth: usize) -> Spectrum {
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
        let photons_per_iteration = params.find_one_int("photonsperiteration", 0) as usize;
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

/// Final image pixel.
#[derive(Default)]
struct SPPMPixel {
    /// Search radius for photons.
    radius: Float,
    /// Estimated average radiance visible over the extent of the pixel including time when
    /// shutter is open to account for depth-of-field.
    ld: Spectrum,
    /// Geometric informtion for the visible point.
    vp: VisiblePoint,
    /// The sum of the product of BSDF values with particle weights
    phi: [AtomicFloat; SPECTRUM_SAMPLES],
    /// Number of photons that have contributed during current interation.
    m: AtomicI32,
    /// Number of photons that have contributed to the point after i^th iteration.
    n: Float,
    /// The sum of products of photons with BSDF values.
    tau: Spectrum,
}

impl SPPMPixel {
    /// Create a new `SPPMPixel` with initial search radius and default field values.
    ///
    /// * `radius` - The search radius for photons.
    fn new(radius: Float) -> Self {
        Self {
            radius,
            ..Default::default()
        }
    }
}

#[derive(Default, Clone)]
struct VisiblePoint {
    /// Point along the camera path to look for nearby photons.
    p: Point3f,
    /// Outgoing direction used to compute reflected radiance.
    wo: Vector3f,
    /// The BSDF at the point used to compute reflected radiance.
    bsdf: Option<BSDF>,
    /// Throughput weight used to compute reflected radiance.
    beta: Spectrum,
}

impl VisiblePoint {
    /// Create a new `VisiblePoint`.
    ///
    /// `p`    - Point along the camera path to look for nearby photons.
    /// `wo`   - Outgoing direction used to compute reflected radiance.
    /// `bsdf` - The BSDF at the point used to compute reflected radiance.
    /// `beta` - Throughput weight used to compute reflected radiance.
    fn new(p: Point3f, wo: Vector3f, bsdf: Option<BSDF>, beta: Spectrum) -> Self {
        Self { p, wo, bsdf, beta }
    }
}

/// Returns the coordinates of the voxel in the grid that the given point lies in. The Boolean return value indicates
/// whether the point is inside the grid’s bounds; if it isn’t, the returned coordinates are clamped to be inside the
/// range of valid coordinates.
///
/// * `p`        - Voxel point.
/// * `bounds`   - Grid bounds.
/// * `grid_res` - Grid resolution.
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

/// Hashes the coordinates of the voxel, returning an index into the grid array.
///
/// * `p`         - Voxel point.
/// * `hash_size` - Hash size (the total number of pixels).
#[inline(always)]
fn hash(p: &Point3i, hash_size: usize) -> usize {
    let (x, _) = p.x.overflowing_mul(73856093);
    let (y, _) = p.y.overflowing_mul(19349663);
    let (z, _) = p.z.overflowing_mul(83492791);

    (x ^ y ^ z) as usize % hash_size
}

/// Compute tile bounds for SPPM tile.
fn get_tile_bounds(tile_x: usize, tile_y: usize, tile_size: usize, pixel_bounds: &Bounds2i) -> Bounds2i {
    let x0 = pixel_bounds.p_min.x + (tile_x * tile_size) as Int;
    let x1 = min(x0 + tile_size as Int, pixel_bounds.p_max.x);
    let y0 = pixel_bounds.p_min.y + (tile_y * tile_size) as Int;
    let y1 = min(y0 + tile_size as Int, pixel_bounds.p_max.y);
    Bounds2i::new(Point2i::new(x0, y0), Point2i::new(x1, y1))
}

/// Generate SPPM visible points.
fn generate_visible_points(
    pixels: &mut [SPPMPixel],
    iter: usize,
    n_tiles: &Point2<usize>,
    tile_size: usize,
    tile_count: usize,
    scene: &Scene,
    sampler: &HaltonSampler,
    pixel_bounds: &Bounds2i,
    inv_sqrt_spp: Float,
    max_depth: usize,
    camera: ArcCamera,
) {
    let n_threads = OPTIONS.threads();
    let camera = &camera;

    thread::scope(|scope| {
        let (tx_collector, rx_collector) =
            crossbeam_channel::bounded::<(usize, Spectrum, Option<VisiblePoint>)>(n_threads);
        let (tx_worker, rx_worker) = crossbeam_channel::bounded::<usize>(n_threads);

        // Spawn collector thread.
        scope.spawn(move || {
            for (pixel_offset, ld, vp) in rx_collector.iter() {
                pixels[pixel_offset].ld += ld;
                if let Some(vp) = vp {
                    pixels[pixel_offset].vp = vp;
                }
            }
        });

        // Spawn worker threads.
        for _ in 0..n_threads {
            let rx_worker = rx_worker.clone();
            let tx_collector = tx_collector.clone();
            scope.spawn(move || {
                for tile_idx in rx_worker.iter() {
                    // Follow camera paths for `tile_`in image for SPPM.
                    let tile_x = tile_idx % n_tiles.x;
                    let tile_y = tile_idx / n_tiles.x;

                    // Get sampler instance for tile.
                    let mut tile_sampler = Sampler::clone(sampler, tile_idx as u64);

                    // Compute `tile_bounds` for SPPM tile.
                    let tile_bounds = get_tile_bounds(tile_x, tile_y, tile_size, &pixel_bounds);

                    for p_pixel in tile_bounds {
                        let p_pixel_o = Point2i::from(p_pixel - pixel_bounds.p_min);
                        let pixel_offset =
                            (p_pixel_o.x + p_pixel_o.y * (pixel_bounds.p_max.x - pixel_bounds.p_min.x)) as usize;

                        if let Some((ld, vp)) = generate_visible_point(
                            iter,
                            &p_pixel,
                            scene,
                            &mut tile_sampler,
                            inv_sqrt_spp,
                            max_depth,
                            Arc::clone(camera),
                        ) {
                            tx_collector.send((pixel_offset, ld, vp)).unwrap();
                        }
                    }
                }
            });
        }
        drop(rx_worker); // Drop extra since we've cloned one for each woker.
        drop(tx_collector);

        // Send work.
        for tile_idx in 0..tile_count {
            tx_worker.send(tile_idx).unwrap();
        }
    });
}

/// Generate single SPPM visible point.
fn generate_visible_point(
    iter: usize,
    p_pixel: &Point2i,
    scene: &Scene,
    tile_sampler: &mut ArcSampler,
    inv_sqrt_spp: Float,
    max_depth: usize,
    camera: ArcCamera,
) -> Option<(Spectrum, Option<VisiblePoint>)> {
    let camera_sample = {
        // Prepare `tile_sampler` for `p_pixel`.
        let tile_sampler = Arc::get_mut(tile_sampler).unwrap();
        tile_sampler.start_pixel(&p_pixel);
        tile_sampler.set_sample_number(iter);

        // Generate camera sample for `p_pixel`.
        tile_sampler.get_camera_sample(&p_pixel)
    };

    // Generate camera ray for pixel for SPPM.
    let (mut ray, beta) = camera.generate_ray_differential(&camera_sample);
    if beta == 0.0 {
        return None;
    }
    let mut beta = Spectrum::new(beta);
    ray.scale_differentials(inv_sqrt_spp);

    // Follow camera ray path until a visible point is created.
    let mut pixel_ld = Spectrum::ZERO;
    let mut pixel_vp: Option<VisiblePoint> = None;

    // Get `SPPMPixel` for `p_pixel`.
    let mut specular_bounce = false;
    let mut depth = 0;
    while depth < max_depth {
        let isect = scene.intersect(&mut ray);
        if isect.is_none() {
            // Accumulate light contributions for ray with no intersection.
            for light in scene.lights.iter() {
                pixel_ld += beta * light.le(&ray);
            }
            break;
        }
        let mut isect = isect.unwrap();

        // Process SPPM camera ray intersection.

        // Compute BSDF at SPPM camera ray intersection.
        let mut bsdf: Option<BSDF> = None;
        let mut bssrdf: Option<BSDF> = None;
        isect.compute_scattering_functions(&mut ray, true, TransportMode::Radiance, &mut bsdf, &mut bssrdf);
        if bsdf.is_none() {
            ray = isect.spawn_ray(&ray.d);
            continue; // No need to decrement depth.
        }
        let bsdf = bsdf.unwrap();

        // Accumulate direct illumination at SPPM camera ray intersection.
        let wo = -ray.d;
        if depth == 0 || specular_bounce {
            pixel_ld += beta * isect.le(&wo);
        }

        let it = Interaction::Surface { si: isect };
        pixel_ld += beta * uniform_sample_one_light(&it, Some(&bsdf), scene, tile_sampler, false, None);

        // Need to move the `isect` back out of `it`.
        let isect = match it {
            Interaction::Surface { si } => si,
            _ => unreachable!(), // isect was SurfaceInteraction. So this should not be possible.
        };

        // Possibly create visible point and end camera path.
        let is_diffuse =
            bsdf.num_components(BxDFType::BSDF_DIFFUSE | BxDFType::BSDF_REFLECTION | BxDFType::BSDF_TRANSMISSION) > 0;
        let is_glossy =
            bsdf.num_components(BxDFType::BSDF_GLOSSY | BxDFType::BSDF_REFLECTION | BxDFType::BSDF_TRANSMISSION) > 0;

        let bsdf = if is_diffuse || (is_glossy && depth == max_depth - 1) {
            pixel_vp = Some(VisiblePoint::new(isect.hit.p, wo, Some(bsdf), beta));
            pixel_vp.as_ref().unwrap().bsdf.as_ref().unwrap() // just so we can use bsdf later
        } else {
            &bsdf
        };

        // Spawn ray from SPPM camera path vertex.
        if depth < max_depth - 1 {
            let tile_sampler = Arc::get_mut(tile_sampler).unwrap();
            let BxDFSample { f, pdf, wi, bxdf_type } = bsdf.sample_f(&wo, &tile_sampler.get_2d(), BxDFType::all());
            if pdf == 0.0 || f.is_black() {
                break;
            }
            specular_bounce = (bxdf_type & BxDFType::BSDF_SPECULAR) > BxDFType::BSDF_NONE;
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

    Some((pixel_ld, pixel_vp))
}

/// Add visible points to SPPM grid.
fn add_visible_points_to_grid<'p>(
    pixels: &'p [SPPMPixel],
    grid_bounds: &Bounds3f,
    grid_res: &[Int; 3],
    hash_size: usize,
    grid: &[Arc<Mutex<VecDeque<&'p SPPMPixel>>>],
) {
    let n_threads = OPTIONS.threads();

    thread::scope(|scope| {
        let (tx_worker, rx_worker) = crossbeam_channel::bounded::<(usize, &SPPMPixel)>(4096);

        // Spawn worker threads.
        for _ in 0..n_threads {
            let rx_worker = rx_worker.clone();
            scope.spawn(move || {
                for (index, pixel) in rx_worker.iter() {
                    if !pixel.vp.beta.is_black() {
                        // Add pixel's visible point to applicable grid cells.
                        let radius = Vector3f::new(pixel.radius, pixel.radius, pixel.radius);

                        let (_, p_min) = to_grid(&(pixel.vp.p - radius), &grid_bounds, &grid_res);
                        let (_, p_max) = to_grid(&(pixel.vp.p + radius), &grid_bounds, &grid_res);

                        for z in p_min.z..=p_max.z {
                            for y in p_min.y..=p_max.y {
                                for x in p_min.x..=p_max.x {
                                    // Add visible point to grid cell `(x, y, z)`.
                                    let h = hash(&Point3i::new(x, y, z), hash_size);
                                    let mut grid_h = grid[h].lock().unwrap();
                                    (*grid_h).push_front(&pixels[index]);
                                }
                            }
                        }
                    }
                }
            });
        }
        drop(rx_worker); // Drop extra since we've cloned one for each woker.

        // Send work.
        for (index, pixel) in pixels.iter().enumerate() {
            tx_worker.send((index, pixel)).unwrap();
        }
    });
}

/// Compute grid bounds for SPPM visible points.
fn compute_grid_bounds(pixels: &[SPPMPixel]) -> (Bounds3f, Float) {
    let mut grid_bounds = Bounds3f::default();
    let mut max_radius = 0.0;
    for pixel in pixels.iter() {
        if pixel.vp.beta.is_black() {
            continue;
        }
        let vp_bound = Bounds3f::from(pixel.vp.p).expand(pixel.radius);
        grid_bounds = grid_bounds.union(&vp_bound);
        max_radius = max(max_radius, pixel.radius);
    }
    (grid_bounds, max_radius)
}

/// Compute resolution of SPPM grid in each dimension.
fn compute_grid_resolution(grid_bounds: &Bounds3f, max_radius: Float) -> [Int; 3] {
    let diag = grid_bounds.diagonal();
    let max_diag = diag.max_component();
    let base_grid_res = (max_diag / max_radius) as Int;
    assert!(base_grid_res > 0);

    let mut grid_res: [Int; 3] = [0; 3];
    for i in 0..3 {
        grid_res[i] = max((base_grid_res as Float * diag[i] / max_diag) as Int, 1);
    }
    grid_res
}

/// Trace photons and accumulate contributions.
fn trace_photons<'p>(
    iter: usize,
    photons_per_iteration: usize,
    scene: &Scene,
    camera_data: &CameraData,
    light_distr: &Option<Distribution1D>,
    max_depth: usize,
    grid_bounds: &Bounds3f,
    grid_res: &[Int; 3],
    hash_size: usize,
    grid: &[Arc<Mutex<VecDeque<&'p SPPMPixel>>>],
) {
    let n_threads = OPTIONS.threads();

    thread::scope(|scope| {
        let (tx_worker, rx_worker) = crossbeam_channel::bounded::<usize>(8192);

        // Spawn worker threads.
        for _ in 0..n_threads {
            let rx_worker = rx_worker.clone();
            scope.spawn(move || {
                for photon_index in rx_worker.iter() {
                    trace_photon(
                        iter,
                        photon_index,
                        photons_per_iteration,
                        scene,
                        camera_data,
                        &light_distr,
                        max_depth,
                        &grid_bounds,
                        &grid_res,
                        hash_size,
                        grid,
                    );
                }
            });
        }
        drop(rx_worker); // Drop extra since we've cloned one for each woker.

        // Send work.
        for photon_index in 0..photons_per_iteration {
            tx_worker.send(photon_index).unwrap();
        }
    });
}

/// Trace single photon and accumulate contribution.
fn trace_photon<'p>(
    iter: usize,
    photon_index: usize,
    photons_per_iteration: usize,
    scene: &Scene,
    camera_data: &CameraData,
    light_distr: &Option<Distribution1D>,
    max_depth: usize,
    grid_bounds: &Bounds3f,
    grid_res: &[Int; 3],
    hash_size: usize,
    grid: &[Arc<Mutex<VecDeque<&'p SPPMPixel>>>],
) {
    // Follow photon path for `photon_index`.
    let halton_index = (iter * photons_per_iteration + photon_index) as u64;
    let mut halton_dim = 0;

    // Choose light to shoot photon from.
    let light_sample = radical_inverse(halton_dim, halton_index);
    halton_dim += 1;

    let mut beta = Spectrum::ZERO;
    let mut photon_ray = Ray::default();
    if let Some(distr) = light_distr {
        let (light_num, light_pdf, _u_remapped) = distr.sample_discrete(light_sample);
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

        // Generate `photon_ray` from light source and initialize `beta`.
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
    while depth < max_depth {
        let isect = scene.intersect(&mut photon_ray);
        if isect.is_none() {
            break;
        }
        let mut isect = isect.unwrap();

        if depth > 0 {
            // Add photon contribution to nearby visible points.
            let (result, photon_grid_index) = to_grid(&isect.hit.p, &grid_bounds, &grid_res);
            if result {
                let h = hash(&photon_grid_index, hash_size);

                // Add photon contribution to visible points in `grid[h]`.
                let grid_h = grid[h].lock().unwrap();
                for pixel in (*grid_h).iter() {
                    let radius = pixel.radius;
                    if pixel.vp.p.distance_squared(isect.hit.p) > radius * radius {
                        continue;
                    }

                    // Update `pixel`, `Phi` and `M` for nearby photon.
                    if let Some(bsdf) = pixel.vp.bsdf.as_ref() {
                        let wi = -photon_ray.d;
                        let phi = beta * bsdf.f(&pixel.vp.wo, &wi, BxDFType::all());
                        for i in 0..SPECTRUM_SAMPLES {
                            pixel.phi[i].add(phi[i], Ordering::AcqRel);
                        }
                        pixel.m.fetch_add(1, Ordering::AcqRel);
                    }
                }
            }
        }
        // Sample new photon ray direction.

        // Compute BSDF at photon intersection point.
        let mut bsdf: Option<BSDF> = None;
        let mut bssrdf: Option<BSDF> = None;
        isect.compute_scattering_functions(&mut photon_ray, true, TransportMode::Importance, &mut bsdf, &mut bssrdf);
        if bsdf.is_none() {
            photon_ray = isect.spawn_ray(&photon_ray.d);
            continue; // No need to decrement depth.
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

/// Update pixel values from this pass' photons.
fn update_pixels(pixels: &mut Vec<SPPMPixel>) {
    let n_threads = OPTIONS.threads();

    thread::scope(|scope| {
        let (tx_worker, rx_worker) = crossbeam_channel::bounded::<&mut SPPMPixel>(4096);

        // Spawn worker threads.
        for _ in 0..n_threads {
            let rx_worker = rx_worker.clone();
            scope.spawn(move || {
                for p in rx_worker.iter() {
                    let m = p.m.load(Ordering::Acquire);
                    if m > 0 {
                        // Update pixel photon count, search radius, and `tau` from photons.
                        const GAMMA: Float = 2.0 / 3.0;
                        let n_new = p.n + GAMMA * m as Float;
                        let r_new = p.radius * (n_new / (p.n + m as Float)).sqrt();

                        let mut phi = Spectrum::ZERO;
                        for j in 0..SPECTRUM_SAMPLES {
                            phi[j] = p.phi[j].load(Ordering::Acquire);
                        }

                        p.tau = (p.tau + p.vp.beta * phi) * (r_new * r_new) / (p.radius * p.radius);
                        p.n = n_new;
                        p.radius = r_new;
                        p.m.store(0, Ordering::Release);

                        for j in 0..SPECTRUM_SAMPLES {
                            p.phi[j].store(0.0, Ordering::Release);
                        }
                    }
                    // Reset `VisiblePoint` in pixel.
                    p.vp.beta = Spectrum::ZERO;
                    p.vp.bsdf = None;
                }
            });
        }
        drop(rx_worker); // Drop extra since we've cloned one for each woker.

        // Send work.
        for pixel in pixels.iter_mut() {
            tx_worker.send(pixel).unwrap();
        }
    });
}

/// Store SPPM image in film and write image.
fn write_sppm_image(
    pixels: &[SPPMPixel],
    iter: usize,
    pixel_bounds: &Bounds2i,
    photons_per_iteration: usize,
    camera_data: &CameraData,
    elapsed_progress: Duration,
) {
    let n_pixels = pixels.len();

    let x0 = pixel_bounds.p_min.x;
    let x1 = pixel_bounds.p_max.x;
    let np = (iter + 1) as u64 * photons_per_iteration as u64;

    let mut image = vec![Spectrum::ZERO; n_pixels];
    let mut offset = 0;
    for y in pixel_bounds.p_min.y..pixel_bounds.p_max.y {
        for x in x0..x1 {
            // Compute radiance `L` for SPPM pixel `pixel`.
            let pixel_index = (y - pixel_bounds.p_min.y) * (x1 - x0) + (x - x0);
            let pixel = &pixels[pixel_index as usize];
            let mut l = pixel.ld / (iter + 1) as Float;
            l += pixel.tau / (np as Float * PI * pixel.radius * pixel.radius);
            image[offset] = l;
            offset += 1;
        }
    }
    camera_data.film.set_image(&image);
    camera_data.film.write_image(1.0);

    // Write SPPM radius image, if requested.
    if let Some(sppm_radius) = OPTIONS.sppm_radius.as_ref() {
        let mut rimg: Vec<Float> = vec![0.0; 3 * n_pixels];
        let mut minrad: Float = 1e30;
        let mut maxrad: Float = 0.0;
        for y in pixel_bounds.p_min.y..pixel_bounds.p_max.y {
            for x in x0..x1 {
                let pixel_index = (y - pixel_bounds.p_min.y) * (x1 - x0) + (x - x0);
                let p = &pixels[pixel_index as usize];
                minrad = min(minrad, p.radius);
                maxrad = max(maxrad, p.radius);
            }
        }
        debug!(
            "Iterations: {} ({:2} s) radius range: {} - {}",
            iter + 1,
            elapsed_progress.as_secs_f32(),
            minrad,
            maxrad
        );
        let mut offset = 0;
        for y in pixel_bounds.p_min.y..pixel_bounds.p_max.y {
            for x in x0..x1 {
                let pixel_index = (y - pixel_bounds.p_min.y) * (x1 - x0) + (x - x0);
                let p = &pixels[pixel_index as usize];
                let v = 1.0 - (p.radius - minrad) / (maxrad - minrad);
                rimg[offset] = v;
                rimg[offset + 1] = v;
                rimg[offset + 2] = v;
                offset += 3;
            }
        }
        if let Err(err) = write_image(sppm_radius, &rimg, &pixel_bounds) {
            error!("Error writing output image sppm_radius.png. {:}.", err);
        }
    }
}
