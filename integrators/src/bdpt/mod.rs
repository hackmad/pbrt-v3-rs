//! Bi-directional Path Tracer

#![allow(dead_code)]

use core::app::OPTIONS;
use core::camera::*;
use core::film::Film;
use core::film::FilmTile;
use core::geometry::*;
use core::integrator::*;
use core::interaction::*;
use core::light::*;
use core::light_distrib::*;
use core::material::*;
use core::paramset::*;
use core::pbrt::*;
use core::reflection::*;
use core::sampler::*;
use core::sampling::*;
use core::scene::*;
use core::spectrum::*;
use filters::*;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::thread;

mod bdpt_interaction;
mod common;
mod vertex;

use common::*;
use vertex::*;

/// Implements bi-directional path tracing integrator.
pub struct BDPTIntegrator {
    /// The sampler.
    sampler: ArcSampler,

    /// The camerea.
    camera: ArcCamera,

    /// Maximum path depth.
    max_depth: usize,

    /// Used to visual the strategies for debug purposes.
    visualize_strategies: bool,

    /// Used to visual the weights for debug purposes.
    visualize_weights: bool,

    /// Pixel bounds for the image.
    pixel_bounds: Bounds2i,

    /// Light sampling strategy.
    light_sample_strategy: LightSampleStategy,

    /// Stores 1 / samples per pixel.
    inv_sample_count: Float,
}

impl BDPTIntegrator {
    /// Create a new `BDPTIntegrator`.
    ///
    /// * `sampler`               - The sampler.
    /// * `camera`                - The camera.
    /// * `max_depth`             - The maximum depth.
    /// * `visualize_strategies`  - Used to visual the strategies for debug purposes.
    /// * `visualize_weights`     - Used to visual the weights for debug purposes.
    /// * `pixel_bounds`          - Pixel bounds for the image.
    /// * `light_sample_strategy` - Light sampling strategy (default to Power)
    pub fn new(
        sampler: ArcSampler,
        camera: ArcCamera,
        max_depth: usize,
        visualize_strategies: bool,
        visualize_weights: bool,
        pixel_bounds: Bounds2i,
        light_sample_strategy: LightSampleStategy,
    ) -> Self {
        let inv_sample_count = 1.0 / sampler.get_data().samples_per_pixel as Float;

        Self {
            sampler,
            camera,
            max_depth,
            visualize_strategies,
            visualize_weights,
            pixel_bounds,
            light_sample_strategy,
            inv_sample_count,
        }
    }

    /// Render an image tile.
    ///
    /// * `tile_idx`       - Unique tile index.
    /// * `n_tiles`        - Number of tiles in (x, y) direction.
    /// * `scene`          - Scene.
    /// * `sample_bounds`  - Sample bounds.
    /// * `light_distr`    - Light distribution.
    /// * `light_to_index` - Map of light distribution indices.
    /// * `weight_films`   - Used for debug visualization.
    /// * `film`           - The film.
    fn render_tile(
        &self,
        tile_idx: usize,
        n_tiles: Point2<usize>,
        scene: &Scene,
        sample_bounds: Bounds2i,
        light_distr: ArcLightDistribution,
        light_to_index: Arc<HashMap<usize, usize>>,
        weight_films: Arc<Mutex<Vec<Option<Film>>>>,
        film: &Film,
    ) -> FilmTile {
        // Get the x and y tile indices.
        let tile_x = tile_idx % n_tiles.x;
        let tile_y = tile_idx / n_tiles.x;

        // Get camera data.
        let camera_data = self.camera.get_data();

        // Get sampler instance for tile.
        let mut tile_sampler = Sampler::clone(&*self.sampler, tile_idx as u64);

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

            if !self.pixel_bounds.contains_exclusive(&pixel) {
                continue;
            }

            loop {
                // Generate a single sample using BDPT.
                let p_film = {
                    let sampler = Arc::get_mut(&mut tile_sampler).unwrap();
                    Point2f::from(pixel) + sampler.get_2d()
                };

                // Trace the camera subpath.
                let mut camera_vertices: Vec<Vertex> = vec![Vertex::default(); self.max_depth + 2];
                let mut light_vertices: Vec<Vertex> = vec![Vertex::default(); self.max_depth + 1];
                let n_camera = generate_camera_subpath(
                    scene,
                    &mut tile_sampler,
                    self.max_depth + 2,
                    &self.camera,
                    &p_film,
                    &mut camera_vertices,
                );

                // Get a distribution for sampling the light at the start of the light subpath. Because
                // the light path follows multiple bounces, basing the sampling distribution on any of
                // the vertices of the camera path is unlikely to be a good strategy. We use the
                // PowerLightDistribution by default here, which doesn't use the point passed to it.
                let light_distr = light_distr.lookup(&camera_vertices[0].it.p()).unwrap();

                // Now trace the light subpath
                let n_light = generate_light_subpath(
                    scene,
                    &mut tile_sampler,
                    self.max_depth + 1,
                    camera_vertices[0].it.time(),
                    Arc::clone(&light_distr),
                    Arc::clone(&light_to_index),
                    &mut light_vertices,
                );

                // Execute all BDPT connection strategies.
                let mut l = Spectrum::ZERO;
                for t in 1..=n_camera {
                    for s in 0..=n_light {
                        let depth = t as isize + s as isize - 2;
                        if (s == 1 && t == 1) || depth < 0 || depth > self.max_depth as isize {
                            continue;
                        }

                        // Execute the $(s, t)$ connection strategy and update `L`.
                        let (l_path, p_film_new, mis_wt) = connect_bdpt(
                            scene,
                            &mut light_vertices,
                            &mut camera_vertices,
                            s,
                            t,
                            Arc::clone(&light_distr),
                            Arc::clone(&light_to_index),
                            &self.camera,
                            &mut tile_sampler,
                        );

                        info!("Connect bdpt s: {s}, t: {t}, l_path: {l_path}, mis_weight: {mis_wt}");
                        if self.visualize_strategies || self.visualize_weights {
                            let mut value = Spectrum::ZERO;
                            if self.visualize_strategies {
                                if mis_wt != 0.0 {
                                    value = l_path / mis_wt;
                                }
                            }
                            if self.visualize_weights {
                                value = l_path;
                            }
                            let weight_films = weight_films.lock().unwrap();
                            if let Some(wf) = &(*weight_films)[buffer_index(s, t)] {
                                wf.add_splat(&p_film_new.unwrap_or(p_film), &value);
                            }
                        }
                        if t != 1 {
                            l += l_path;
                        } else {
                            film.add_splat(&p_film_new.unwrap_or(p_film), &l_path);
                        }
                    }
                }

                info!("Add film sample pFilm: {p_film}, L: {l}, (y: {})", l.y());
                film_tile.add_sample(p_film, l, 1.0);

                {
                    let sampler = Arc::get_mut(&mut tile_sampler).unwrap();
                    if !sampler.start_next_sample() {
                        break;
                    }
                }
            }
        }

        film_tile
    }
}

impl Integrator for BDPTIntegrator {
    /// Preprocess the scene.
    ///
    /// * `scene` - The scene
    fn preprocess(&mut self, _scene: &Scene) {}

    /// Render the scene.
    ///
    /// * `scene` - The scene.
    fn render(&self, scene: &Scene) {
        let light_distribution = create_light_sample_distribution(self.light_sample_strategy, scene);

        // Compute a reverse mapping from light pointers to offsets into the scene lights vector (and,
        // equivalently, offsets into light_distr). Added after book text was finalized; this is critical
        // to reasonable performance with 100s+ of light sources.
        let mut light_to_index = HashMap::new();
        for (i, light) in scene.lights.iter().enumerate() {
            light_to_index.insert(light.get_id(), i);
        }
        let light_to_index: Arc<HashMap<usize, usize>> = Arc::new(light_to_index);

        // Partition the image into tiles.
        let camera_data = self.camera.get_data();
        let sample_bounds = camera_data.film.get_sample_bounds();
        let sample_extent = sample_bounds.diagonal();
        let tile_size = OPTIONS.tile_size;
        let n_tiles = Point2::new(
            ((sample_extent.x + tile_size as Int - 1) / tile_size as Int) as usize,
            ((sample_extent.y + tile_size as Int - 1) / tile_size as Int) as usize,
        );
        let tile_count = n_tiles.x * n_tiles.y;

        // Allocate buffers for debug visualization.
        let buffer_count = (1 + self.max_depth) * (6 + self.max_depth) / 2;
        let mut weight_films: Vec<Option<Film>> = (0..buffer_count).map(|_| None).collect();

        let mut weight_film_count = 0;
        if self.visualize_strategies || self.visualize_weights {
            for depth in 0..=self.max_depth {
                for s in 0..=depth + 2 {
                    let t = depth + 2 - s;
                    if t == 0 || (s == 1 && t == 1) {
                        continue;
                    }

                    let filename = format!("bdpt_d{depth:02}_s{s:02}_t{t:02}.exr");

                    weight_films[buffer_index(s, t)] = Some(Film::new(
                        &camera_data.film.full_resolution,
                        &Bounds2f::new(Point2f::new(0.0, 0.0), Point2f::new(1.0, 1.0)),
                        Arc::new(BoxFilter::from(&ParamSet::new())),
                        camera_data.film.diagonal * 1000.0,
                        &filename,
                        Some(1.0),
                        None,
                    ));
                    weight_film_count += 1;
                }
            }
        }
        let weight_films = Arc::new(Mutex::new(weight_films));

        let progress = create_progress_bar(tile_count as u64 + weight_film_count as u64 + 1);
        progress.set_message("Rendering scene");

        // Render and write the output image to disk.
        if !scene.lights.is_empty() {
            let n_threads = OPTIONS.threads();

            thread::scope(|scope| {
                let (tx, rx) = crossbeam_channel::bounded(n_threads);

                // Spawn worker threads.
                for _ in 0..n_threads {
                    let rxc = rx.clone();
                    let progress = &progress;
                    let light_distr = &light_distribution;
                    let light_to_index = &light_to_index;
                    let weight_films = &weight_films;
                    scope.spawn(move || {
                        for tile_idx in rxc.iter() {
                            // Render section of image corresponding to `tile`.
                            let film_tile = self.render_tile(
                                tile_idx,
                                n_tiles,
                                scene,
                                sample_bounds,
                                Arc::clone(light_distr),
                                Arc::clone(light_to_index),
                                Arc::clone(weight_films),
                                &camera_data.film,
                            );

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
            });
        } else {
            progress.inc(tile_count as u64);
        }

        // Save final image after rendering.
        progress.set_message("Writing image");
        camera_data.film.write_image(self.inv_sample_count);
        progress.inc(1);

        // Write buffers for debug visualization.
        if self.visualize_strategies || self.visualize_weights {
            progress.set_message("Writing debug images");
            let weight_films = weight_films.lock().unwrap();
            for weight_film in weight_films.iter() {
                if let Some(wf) = weight_film {
                    wf.write_image(self.inv_sample_count);
                    progress.inc(1);
                }
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

impl From<(&ParamSet, ArcSampler, ArcCamera)> for BDPTIntegrator {
    /// Create a `PathIntegrator` from given parameter set and camera.
    ///
    /// * `p` - A tuple containing parameter set and camera.
    fn from(p: (&ParamSet, ArcSampler, ArcCamera)) -> Self {
        let (params, sampler, camera) = p;

        let mut max_depth = params.find_one_int("maxdepth", 5) as usize;
        let visualize_strategies = params.find_one_bool("visualizestrategies", false);
        let visualize_weights = params.find_one_bool("visualizeweights", false);
        if (visualize_strategies || visualize_weights) && max_depth > 5 {
            warn!("visualizestrategies/visualizeweights was enabled, limiting maxdepth to 5");
            max_depth = 5;
        }

        let pb = params.find_int("pixelbounds");
        let np = pb.len();

        let mut pixel_bounds = camera.get_data().film.get_sample_bounds();
        if np > 0 {
            if np != 4 {
                error!("Expected 4 values for 'pixel_bounds' parameter. Got {np}");
            } else {
                pixel_bounds =
                    pixel_bounds.intersect(&Bounds2i::new(Point2i::new(pb[0], pb[1]), Point2i::new(pb[2], pb[3])));
                if pixel_bounds.area() == 0 {
                    error!("Degenerate 'pixel_bounds' specified.");
                }
            }
        }

        let lss = params.find_one_string("lightsamplestrategy", "power".to_owned());
        let light_sample_strategy = LightSampleStategy::from(lss.as_ref());

        Self::new(
            sampler,
            camera,
            max_depth,
            visualize_strategies,
            visualize_weights,
            pixel_bounds,
            light_sample_strategy,
        )
    }
}

/// Generate the camera subpath. It returns the number of vertices in the subpath.
///
/// * `scene`     - The scene.
/// * `sampler`   - The sampler.
/// * `max_depth` - The maximum path depth.
/// * `p_film`    - The point on camera film.
/// * `path`      - The path.
fn generate_camera_subpath<'scene>(
    scene: &'scene Scene,
    sampler: &mut ArcSampler,
    max_depth: usize,
    camera: &ArcCamera,
    p_film: &Point2f,
    path: &mut [Vertex<'scene>],
) -> usize {
    if max_depth == 0 {
        return 0;
    }

    // Sample initial ray for camera subpath.
    let (camera_sample, samples_per_pixel) = {
        let sampler = Arc::get_mut(sampler).unwrap();
        let time = sampler.get_1d();
        let p_lens = sampler.get_2d();
        (
            CameraSample::new(*p_film, p_lens, time),
            sampler.get_data().samples_per_pixel,
        )
    };

    let (mut ray, beta) = camera.generate_ray_differential(&camera_sample);
    let beta = Spectrum::new(beta);
    ray.scale_differentials(1.0 / (samples_per_pixel as Float).sqrt());

    // Generate first vertex on camera subpath and start random walk.
    path[0] = Vertex::create_camera_from_ray(Arc::clone(camera), &ray, beta);

    let PDFResult {
        pos: pdf_pos,
        dir: pdf_dir,
    } = camera.pdf_we(&ray);
    info!("Starting camera subpath. Ray: {ray}, beta {beta}, pdf_pos {pdf_pos}, pdf_dir {pdf_dir}");

    let n_vertices = random_walk(
        scene,
        &mut ray,
        sampler,
        beta,
        pdf_dir,
        max_depth - 1,
        TransportMode::Radiance,
        path,
    );

    n_vertices + 1
}

/// Generate the light subpath. It returns the number of vertices in the subpath.
///
/// * `scene`                - The scene.
/// * `sampler`              - The sampler.
/// * `max_depth`            - The maximum depth.
/// * `time`                 - The time.
/// * `light_distr`          - Light probabilities.
/// * `light_to_distr_index` - Map of light distribution indices.
/// * `path`                 - The path.
fn generate_light_subpath<'scene>(
    scene: &'scene Scene,
    sampler: &mut ArcSampler,
    max_depth: usize,
    time: Float,
    light_distr: Arc<Distribution1D>,
    light_to_index: Arc<HashMap<usize, usize>>,
    path: &mut [Vertex<'scene>],
) -> usize {
    if max_depth == 0 {
        return 0;
    }

    // Sample initial ray for light subpath.
    let (light, light_pdf, sample) = {
        let sampler = Arc::get_mut(sampler).unwrap();
        let (light_num, light_pdf, _u_remapped) = light_distr.sample_discrete(sampler.get_1d());
        let light = &scene.lights[light_num];
        let sample = light.sample_le(&sampler.get_2d(), &sampler.get_2d(), time);
        (light, light_pdf, sample)
    };

    let Le {
        mut ray,
        n_light,
        pdf_pos,
        pdf_dir,
        value: le,
    } = sample;

    if pdf_pos == 0.0 || pdf_dir == 0.0 || le.is_black() {
        return 0;
    }

    // Generate first vertex on light subpath and start random walk.
    path[0] = Vertex::create_light_from_ray_normal(Arc::clone(light), &ray, n_light, le, pdf_pos * light_pdf);
    let beta = le * n_light.abs_dot(&ray.d) / (light_pdf * pdf_pos * pdf_dir);

    info!("Starting light subpath. Ray: {ray}, le {le}, beta {beta}, pdf_pos {pdf_pos}, pdf_dir {pdf_dir}");
    let n_vertices = random_walk(
        scene,
        &mut ray,
        sampler,
        beta,
        pdf_dir,
        max_depth - 1,
        TransportMode::Importance,
        path,
    );

    // Correct subpath sampling densities for infinite area lights.
    if path[0].is_infinite_light() {
        // Set spatial density of `path[1]` for infinite area light.
        if n_vertices > 0 {
            path[1].pdf_fwd = pdf_pos;
            if path[1].is_on_surface() {
                path[1].pdf_fwd *= ray.d.abs_dot(&path[1].it.ng());
            }
        }

        // Set spatial density of `path[0]` for infinite area light.
        path[0].pdf_fwd = infinite_light_density(scene, light_distr, light_to_index, &ray.d);
    }

    n_vertices + 1
}

/// Generates the vertices in a camera/light subpath starting at index 1. It returns the number of
/// vertices generated (this does not include `path[0]` which is generated by either
/// `generate_light_subpath()` or `generate_camera_subpath()`).
///
/// * `scene`     - The scene.
/// * `ray`       - The ray containing the position and outgoing direction previously sampled.
/// * `sampler`   - The sampler.
/// * `beta`      - Path throughtput weight.
/// * `pdf`       - Probability of sampling the ray per unit solid angle of `ray.d`.
/// * `max_depth` - The maximum depth.
/// * `mode`      - The light transport mode.
/// * `path`      - The path.
fn random_walk<'scene>(
    scene: &'scene Scene,
    ray: &mut Ray,
    sampler: &mut ArcSampler,
    beta: Spectrum,
    pdf: Float,
    max_depth: usize,
    mode: TransportMode,
    path: &mut [Vertex<'scene>],
) -> usize {
    if max_depth == 0 {
        return 0;
    };

    let mut bounces = 0;
    // Declare variables for forward and reverse probability densities.
    let mut pdf_fwd = pdf;
    let mut pdf_rev = 0.0;

    loop {
        // Attempt to create the next subpath vertex in `path`.
        info!("Random walk. Bounces {bounces}, beta {beta}, pdf_fwd {pdf_fwd}, pdfRev {pdf_rev}");

        // Trace a ray and sample the medium, if any.
        let isect = scene.intersect(ray);

        let (mut beta, mi) = if let Some(medium) = ray.medium.as_ref() {
            let (sample, mi) = medium.sample(&ray, sampler);
            (beta * sample, mi)
        } else {
            (beta, None)
        };

        if beta.is_black() {
            break;
        }

        let vertex = bounces + 1; // Skip path[0]
        let prev = bounces;
        if let Some(mi) = mi {
            // Record medium interaction in _path_ and compute forward density.
            path[vertex] = Vertex::create_medium(mi, beta, pdf_fwd, &path[prev]);
            bounces += 1;
            if bounces >= max_depth {
                break;
            }

            // Sample direction and compute reverse density at preceding vertex.
            let mi = match &path[vertex].it {
                Interaction::Medium { mi } => mi,
                _ => unreachable!(),
            };

            let sampler = Arc::get_mut(sampler).unwrap();
            let (pdf, wi) = mi.phase.sample_p(&-ray.d, &sampler.get_2d());
            pdf_fwd = pdf;
            pdf_rev = pdf;
            *ray = mi.spawn_ray(&wi);
        } else {
            // Handle surface interaction for path generation.
            if isect.is_none() {
                // Capture escaped rays when tracing from the camera.
                if mode == TransportMode::Radiance {
                    let ei = EndpointInteraction::light_from_ray(ray);
                    path[vertex] = Vertex::create_light_from_endpoint_interaction(ei, beta, pdf_fwd);
                    bounces += 1;
                }
                break;
            }
            let mut isect = isect.unwrap();

            // Compute scattering functions for `mode` and skip over medium boundaries.
            let mut bsdf: Option<BSDF> = None;
            let mut bssrdf: Option<BSDF> = None;
            isect.compute_scattering_functions(ray, true, mode, &mut bsdf, &mut bssrdf);
            if bsdf.is_none() {
                *ray = isect.spawn_ray(&ray.d);
                continue;
            }
            let bsdf = bsdf.unwrap();

            // Initialize `vertex` with surface intersection information.
            path[vertex] = Vertex::create_surface(isect, bsdf, beta, pdf_fwd, &path[prev]);
            bounces += 1;
            if bounces >= max_depth {
                break;
            }

            // Borrow these for later.
            let (isect, bsdf) = match &path[vertex].it {
                Interaction::Surface { si } => (si, path[vertex].bsdf.as_ref().unwrap()),
                _ => unreachable!(),
            };

            // Sample BSDF at current vertex and compute reverse probability.
            let wo = isect.hit.wo;

            let sampler = Arc::get_mut(sampler).unwrap();
            let BxDFSample { f, pdf, wi, bxdf_type } = bsdf.sample_f(&wo, &sampler.get_2d(), BxDFType::all());
            pdf_fwd = pdf;

            info!("Random walk sampled dir {wi} f: {f}, pdf_fwd: {pdf_fwd}");
            if f.is_black() || pdf_fwd == 0.0 {
                break;
            }

            beta *= f * wi.abs_dot(&isect.shading.n) / pdf_fwd;
            info!("Random walk beta now {beta}");
            pdf_rev = bsdf.pdf(&wi, &wo, BxDFType::all());
            if bxdf_type & BxDFType::BSDF_SPECULAR > BxDFType::BSDF_NONE {
                path[vertex].delta = true;
                pdf_rev = 0.0;
                pdf_fwd = 0.0;
            }

            beta *= correct_shading_normal(&isect, &wo, &wi, mode);
            info!("Random walk beta after shading normal correction {beta}");
            *ray = isect.spawn_ray(&wi);
        }

        // Compute reverse area density at preceding vertex.
        path[prev].pdf_rev = path[vertex].convert_density(pdf_rev, &path[prev]);
    }

    bounces
}

/// Computes the generalized form of the geometric term `G` for a path between two points.
///
/// * `scene`   - The scene.
/// * `sampler` - The sampler.
/// * `v0`      - First vertex.
/// * `v1`      - Second vertex.
fn g<'scene>(scene: &'scene Scene, sampler: &mut ArcSampler, v0: &Vertex<'scene>, v1: &Vertex<'scene>) -> Spectrum {
    let mut d = v0.it.p() - v1.it.p();

    let mut g1 = 1.0 / d.length_squared();
    d *= g1.sqrt();

    if v0.is_on_surface() {
        g1 *= v0.ns().abs_dot(&d);
    }
    if v1.is_on_surface() {
        g1 *= v1.ns().abs_dot(&d);
    }

    let vis = VisibilityTester::new(v0.it.get_hit().clone(), v1.it.get_hit().clone());
    g1 * vis.tr(scene, sampler)
}

// Define helper function `remap0` that deals with Dirac delta functions. Used to handle special case
// of Dirac delta functions in the path which have a continuous density of 0 which ar emapped to 1.
//
// * `f` - The value to remap.
fn remap0(f: Float) -> Float {
    if f != 0.0 {
        f
    } else {
        1.0
    }
}

/// Calculates the multiple importance sampling weights. It attempts to produce a complete path
/// by iterating over all alternative strategies that could hypothetically have generated the same
/// input path but with an earlier or later crossover point between the light and camera subpaths.
///
/// * `scene`           - The scene.
/// * `light_vertices`  - The vertices in the light subpath.
/// * `camera_vertices` - The vertices in the camera subpath.
/// * `sampled`         - The sampled vertex.
/// * `s`               - Camera subpath prefixes used by a successful BDPT connection.
/// * `t`               - Light subpath prefixes used by a successful BDPT connection.
/// * `light_pdf`       - Light distribution.
/// * `light_to_index`  - Map of light distribution indices.
#[rustfmt::skip]
fn mis_weight<'scene>(
    scene: &'scene Scene,
    light_vertices: &mut [Vertex<'scene>],
    camera_vertices: &mut [Vertex<'scene>],
    sampled: &Vertex<'scene>,
    s: usize,
    t: usize,
    light_pdf: Arc<Distribution1D>,
    light_to_index: Arc<HashMap<usize, usize>>,
) -> Float {
    if s + t == 2 {
        return 1.0;
    };

    let mut sum_ri = 0.0;

    // Lookup connection vertices and predecessors.
    let qs = s as isize - 1;
    let pt = t as isize - 1;
    let qs_minus = qs - 1;
    let pt_minus = pt - 1;

    // The original implementation uses ScopedAssignment which backs up a value and restores
    // that when the variable goes out of scope. This could be done here by storing a &mut T where
    // T is Clone. However, it also using ScopedAssignement on fields of such T's. Taking a mutable
    // reference to the fields of T will make it so that the original T cannot be accessed until the
    // second one is dropped. So we resort to backing up values ourselves we need and restoring them
    // when we are done.

    // Update sampled vertex for `s=1` and `t=1` strategy.
    let mut backup_qs: Option<Vertex> = None;
    let mut backup_pt: Option<Vertex> = None;
    if s == 1 {
        backup_qs = Some(light_vertices[qs as usize].clone());
        light_vertices[qs as usize] = sampled.clone();
    } else if t == 1 {
        backup_pt = Some(camera_vertices[pt as usize].clone());
        camera_vertices[pt as usize] = sampled.clone();
    }

    // Make connection vertices as non-degenerate.
    let mut backup_pt_delta: Option<bool> = None;
    if pt >= 0 {
        backup_pt_delta = Some(camera_vertices[pt as usize].delta);
        camera_vertices[pt as usize].delta = false;
    }

    let mut backup_qs_delta: Option<bool> = None;
    if qs >= 0 {
        backup_qs_delta = Some(light_vertices [qs as usize].delta);
        light_vertices [qs as usize].delta = false;
    }

    // Update reverse density of vertex `pt_t-1`.
    let mut backup_pt_pdf_rev: Option<Float> = None;
    if pt >= 0 {
        backup_pt_pdf_rev = Some(camera_vertices[pt as usize].pdf_rev);
        camera_vertices[pt as usize].pdf_rev = if s > 0 {
            light_vertices[qs as usize].pdf(
                scene, 
                if qs_minus >= 0 { Some(&light_vertices[qs_minus as usize]) } else { None }, 
                &camera_vertices[pt as usize],
            )
        } else {
            assert!(pt_minus >= 0);
            camera_vertices[pt as usize].pdf_light_origin(
                scene,
                &camera_vertices[pt_minus as usize],
                light_pdf,
                light_to_index,
            )
        };
    }

    // Update reverse density of vertex `pt_t-2`.
    let mut backup_pt_minus_pdf_rev: Option<Float> = None;
    if pt_minus >= 0 {
        backup_pt_minus_pdf_rev = Some(camera_vertices[pt_minus as usize].pdf_rev);
        camera_vertices[pt_minus as usize].pdf_rev = if s > 0 {
            camera_vertices[pt as usize].pdf(
                scene, 
                if qs >= 0 { Some(&light_vertices[qs as usize]) } else { None }, 
                &camera_vertices[pt_minus as usize],
            )
        } else {
            camera_vertices[pt as usize].pdf_light(scene, &camera_vertices[pt_minus as usize])
        };
    }

    // Update reverse density of vertices `qs_s-1` and `qs_s-2`.
    let mut backup_qs_pdf_rev: Option<Float> = None;
    if qs >= 0 {
        backup_qs_pdf_rev = Some(light_vertices[qs as usize].pdf_rev);
        light_vertices[qs as usize].pdf_rev = camera_vertices[pt as usize].pdf(
            scene,
            if pt_minus >= 0 { Some(&camera_vertices[pt_minus as usize]) } else { None },
            &light_vertices[qs as usize],
        );
    }

    let mut backup_qs_minus_pdf_rev: Option<Float> = None;
    if qs_minus >= 0 {
        backup_qs_minus_pdf_rev = Some(light_vertices[qs_minus as usize].pdf_rev);
        light_vertices[qs_minus as usize].pdf_rev = light_vertices[qs as usize].pdf(
            scene,
            if pt >= 0 { Some(&camera_vertices[pt as usize]) } else { None },
            &light_vertices[qs_minus as usize],
        )
    }

    // Consider hypothetical connection strategies along the camera subpath.
    let mut ri = 1.0;
    for i in (1..=t-1).rev() { 
        ri *= remap0(camera_vertices[i].pdf_rev) / remap0(camera_vertices[i].pdf_fwd);

        if !camera_vertices[i].delta && !camera_vertices[i - 1].delta {
            sum_ri += ri;
        }
    }

    // Consider hypothetical connection strategies along the light subpath.
    ri = 1.0;
    for i in (0..s).rev() {
        ri *= remap0(light_vertices[i].pdf_rev) / remap0(light_vertices[i].pdf_fwd);

        let delta_light_vertex = if i > 0 {
            light_vertices[i - 1].delta
        } else {
            light_vertices[0].is_delta_light()
        };

        if !light_vertices[i].delta && !delta_light_vertex {
            sum_ri += ri;
        }
    }

    // Restore snapshots in reverse order of backups.
    if let Some(v) = backup_qs_minus_pdf_rev { light_vertices [qs_minus as usize].pdf_rev = v; }
    if let Some(v) = backup_qs_pdf_rev       { light_vertices [qs       as usize].pdf_rev = v; }
    if let Some(v) = backup_pt_minus_pdf_rev { camera_vertices[pt_minus as usize].pdf_rev = v; }
    if let Some(v) = backup_pt_pdf_rev       { camera_vertices[pt       as usize].pdf_rev = v; }

    if let Some(v) = backup_qs_delta { light_vertices [qs as usize].delta = v; }
    if let Some(v) = backup_pt_delta { camera_vertices[pt as usize].delta = v; }

    if let Some(v) = backup_pt { camera_vertices[pt as usize] = v.clone(); }
    if let Some(v) = backup_qs { light_vertices [qs as usize] = v.clone(); }

    1.0 / (1.0 + sum_ri)
}

/// Attempts to connect 2 subpaths with the given number of vertices and returns the weighted
/// contribution of radiance carried along the resulting subpath.
///
/// * `scene`           - The scene.
/// * `light_vertices`  - The vertices in the light subpath.
/// * `camera_vertices` - The vertices in the camera subpath.
/// * `s`               - Number of vertices s from light subpath.
/// * `t`               - Number of vertices t from camera subpath.
/// * `light_distr`     - Light distribution.
/// * `light_to_index`  - Map of light distribution indices.
/// * `camera`          - The camera.
/// * `sampler`         - The sampler.
fn connect_bdpt<'scene>(
    scene: &'scene Scene,
    light_vertices: &mut [Vertex<'scene>],
    camera_vertices: &mut [Vertex<'scene>],
    s: usize,
    t: usize,
    light_distr: Arc<Distribution1D>,
    light_to_index: Arc<HashMap<usize, usize>>,
    camera: &ArcCamera,
    sampler: &mut ArcSampler,
) -> (Spectrum, Option<Point2f>, Float) {
    let mut l = Spectrum::ZERO;
    let mut p_raster: Option<Point2f> = None;

    // Ignore invalid connections related to infinite area lights.
    if t > 1 && s != 0 && matches!(camera_vertices[t - 1].vertex_type, VertexType::Light) {
        return (l, p_raster, 0.0);
    }

    // Perform connection and write contribution to `L`.
    let mut sampled = Vertex::default();
    if s == 0 {
        // Interpret the camera subpath as a complete path.
        let pt = &camera_vertices[t - 1];
        if pt.is_light() {
            assert!(t >= 2);
            l = pt.le(scene, &camera_vertices[t - 2]) * pt.beta;
        }
        debug_assert!(!l.has_nans());
    } else if t == 1 {
        // Sample a point on the camera and connect it to the light subpath.
        assert!(s >= 1);
        let qs = &light_vertices[s - 1];
        if qs.is_connectible() {
            let SampleResult {
                spectrum: spectrum_wi,
                wi,
                pdf,
                p_raster: pr,
                vis,
            } = {
                let sampler = Arc::get_mut(sampler).unwrap();
                camera.sample_wi(&qs.it, &sampler.get_2d())
            };

            p_raster = pr;

            if pdf > 0.0 && !spectrum_wi.is_black() {
                // Initialize dynamically sampled vertex and _L_ for `t=1` case.
                sampled = Vertex::create_camera_from_hit(Arc::clone(camera), vis.p1.clone(), spectrum_wi / pdf);
                l = qs.beta * qs.f(&sampled, TransportMode::Importance) * sampled.beta;

                if qs.is_on_surface() {
                    l *= wi.abs_dot(qs.ns());
                }
                debug_assert!(!l.has_nans());

                // Only check visibility after we know that the path would make a non-zero contribution.
                if !l.is_black() {
                    l *= vis.tr(scene, sampler);
                }
            }
        }
    } else if s == 1 {
        // Sample a point on a light and connect it to the camera subpath.
        assert!(t >= 1);
        let pt = &camera_vertices[t - 1];
        if pt.is_connectible() {
            let (light, light_pdf, wi, pdf, vis, light_weight) = {
                let sampler = Arc::get_mut(sampler).unwrap();

                let (light_num, light_pdf, _u_remapped) = light_distr.sample_discrete(sampler.get_1d());
                let light = &scene.lights[light_num];

                if let Some(li) = light.sample_li(pt.it.get_hit(), &sampler.get_2d()) {
                    (light, light_pdf, li.wi, li.pdf, Some(li.visibility), li.value)
                } else {
                    (light, light_pdf, Vector3f::ZERO, 0.0, None, Spectrum::ZERO)
                }
            };

            if pdf > 0.0 && !light_weight.is_black() {
                let vis = vis.unwrap();

                let ei = EndpointInteraction::light_from_hit(vis.p1.clone(), Some(Arc::clone(light)));
                sampled = Vertex::create_light_from_endpoint_interaction(ei, light_weight / (pdf * light_pdf), 0.0);

                sampled.pdf_fwd =
                    sampled.pdf_light_origin(scene, pt, Arc::clone(&light_distr), Arc::clone(&light_to_index));

                l = pt.beta * pt.f(&sampled, TransportMode::Radiance) * sampled.beta;
                if pt.is_on_surface() {
                    l *= wi.abs_dot(pt.ns());
                }

                // Only check visibility if the path would carry radiance.
                if !l.is_black() {
                    l *= vis.tr(scene, sampler);
                }
            }
        }
    } else {
        // Handle all other bidirectional connection cases
        assert!(s >= 1 && t >= 1);
        let qs = &light_vertices[s - 1];
        let pt = &camera_vertices[t - 1];
        if qs.is_connectible() && pt.is_connectible() {
            l = qs.beta * qs.f(pt, TransportMode::Importance) * pt.f(qs, TransportMode::Radiance) * pt.beta;

            let qs_f = qs.f(pt, TransportMode::Importance);
            let pt_f = pt.f(qs, TransportMode::Radiance);
            let g1 = g(scene, sampler, qs, pt);
            let d_sq = qs.it.p().distance_squared(pt.it.p());

            info!("General connect s: {s}, t: {t}, qs: {qs}, pt: {pt}, qs.f(pt): {qs_f}, pt.f(qs): {pt_f}, G: {g1}, dist^2: {d_sq}");

            if !l.is_black() {
                l *= g1;
            }
        }
    }

    // Compute MIS weight for connection strategy
    let mis_wt = if l.is_black() {
        0.0
    } else {
        mis_weight(
            scene,
            light_vertices,
            camera_vertices,
            &sampled,
            s,
            t,
            light_distr,
            light_to_index,
        )
    };
    info!("MIS weight for (s,t) = ({s}, {t}) connection: {mis_wt}");
    debug_assert!(!mis_wt.is_nan());
    l *= mis_wt;

    (l, p_raster, mis_wt)
}

#[inline]
fn buffer_index(s: usize, t: usize) -> usize {
    let above = s as isize + t as isize - 2;
    (s as isize + above * (5 + above) / 2) as usize
}
