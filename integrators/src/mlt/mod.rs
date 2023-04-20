//! MLTIntegrator

use crate::bdpt::*;
use core::app::OPTIONS;
use core::camera::*;
use core::geometry::*;
use core::integrator::*;
use core::paramset::ParamSet;
use core::pbrt::*;
use core::report_stats;
use core::rng::*;
use core::sampler::*;
use core::sampling::*;
use core::scene::*;
use core::spectrum::*;
use core::stats::*;
use std::collections::HashMap;
use std::sync::Arc;
use std::thread;

mod mlt_sampler;
use mlt_sampler::*;

/// Implements Multiplexed Metropolis Light Transport with Primary Sample Space Sampler to render images and
/// Bidirectional Path Tracing.
pub struct MLTIntegrator {
    /// The camera.
    camera: ArcCamera,

    /// Maximum path depth.
    max_depth: usize,

    /// Sets the number of bootstrapping samples to use to seed the iterations and compute the integral of the scalar
    /// contribution function.
    n_bootstrap: usize,

    /// Number of parallel Markov chains.
    n_chains: usize,

    /// Number of iterations that MLT (on average!) spends in each pixel.
    mutations_per_pixel: usize,

    /// Controls the size of “small step” mutations.
    sigma: Float,

    /// Probability of taking a “large step” mutation.
    large_step_probability: Float,
}

impl MLTIntegrator {
    /// Create a new `MLTIntegrator`.
    ///
    /// * `camera`                 - The camera.
    /// * `max_depth`              - Maximum path depth.
    /// * `n_bootstrap`            - Sets the number of bootstrapping samples to use to seed the iterations and compute
    ///                              the integral of the scalar contribution function.
    /// * `n_chains`               - Number of parallel Markov chains.
    /// * `mutations_per_pixel`    - Number of iterations that MLT (on average!) spends in each pixel.
    /// * `sigma`                  - Controls the size of “small step” mutations.
    /// * `large_step_probability` - Probability of taking a “large step” mutation.
    pub fn new(
        camera: ArcCamera,
        max_depth: usize,
        n_bootstrap: usize,
        n_chains: usize,
        mutations_per_pixel: usize,
        sigma: Float,
        large_step_probability: Float,
    ) -> Self {
        Self {
            camera,
            max_depth,
            n_bootstrap,
            n_chains,
            mutations_per_pixel,
            sigma,
            large_step_probability,
        }
    }

    /// Computes the radiance `L(X)` for a vector of sample values provided by an `MLTSampler`. It returns the raster
    /// position of path and radiance value.
    ///
    /// * `scene` - The scene.
    /// * `light_distr`          - Light probabilities.
    /// * `light_to_distr_index` - Map of light distribution indices.
    /// * `sampler`              - Primary Sample Space Sampler.
    /// * `depth`                - Path depth.
    fn l(
        &self,
        scene: &Scene,
        light_distr: Arc<Distribution1D>,
        light_to_index: Arc<HashMap<usize, usize>>,
        sampler: &mut MLTSampler,
        depth: usize,
    ) -> (Point2f, Spectrum) {
        sampler.start_stream(CAMERA_STREAM_INDEX);

        // Determine the number of available strategies and pick a specific one.
        let (n_strategies, s, t) = if depth == 0 {
            (1, 0, 2)
        } else {
            let n_strategies = depth + 2;
            let s = ((sampler.get_1d() * n_strategies as Float) as isize).min(n_strategies as isize - 1) as usize;
            (n_strategies, s, n_strategies - s)
        };

        // Generate a camera subpath with exactly `t` vertices.
        let mut camera_vertices: Vec<Vertex> = vec![Vertex::default(); t];
        let sample_bounds = Bounds2f::from(self.camera.get_data().film.get_sample_bounds());
        let p_raster = sample_bounds.lerp(&sampler.get_2d());

        let n_camera = generate_camera_subpath(scene, sampler, t, &self.camera, &p_raster, &mut camera_vertices);
        if n_camera != t {
            return (p_raster, Spectrum::ZERO);
        }

        // Generate a light subpath with exactly `s` vertices.
        sampler.start_stream(LIGHT_STREAM_INDEX);
        let mut light_vertices: Vec<Vertex> = vec![Vertex::default(); s];
        let n_light = generate_light_subpath(
            scene,
            sampler,
            s,
            camera_vertices[0].it.time(),
            Arc::clone(&light_distr),
            Arc::clone(&light_to_index),
            &mut light_vertices,
        );
        if n_light != s {
            return (p_raster, Spectrum::ZERO);
        }

        // Execute connection strategy and return the radiance estimate.
        sampler.start_stream(CONNECTION_STREAM_INDEX);
        let (l_path, p_raster_new, _mis_wt) = connect_bdpt(
            scene,
            &mut light_vertices,
            &mut camera_vertices,
            s,
            t,
            Arc::clone(&light_distr),
            Arc::clone(&light_to_index),
            &self.camera,
            sampler,
        );
        (p_raster_new.unwrap_or(p_raster), l_path * n_strategies as Float)
    }
}

impl Integrator for MLTIntegrator {
    /// Render the scene.
    ///
    /// * `scene` - The scene.
    fn render(&self, scene: &Scene) {
        let n_threads = OPTIONS.threads();

        let light_distr = compute_light_power_distribution(scene).map(Arc::new).unwrap();

        // Compute a reverse mapping from light pointers to offsets into the scene lights vector (and, equivalently,
        // offsets into lightDistr). Added after book text was finalized; this is critical to reasonable performance
        // with 100s+ of light sources.
        let mut light_to_index = HashMap::new();
        for (i, light) in scene.lights.iter().enumerate() {
            light_to_index.insert(light.get_id(), i);
        }
        let light_to_index: Arc<HashMap<usize, usize>> = Arc::new(light_to_index);

        // Generate bootstrap samples and compute normalization constant `b`.
        let n_bootstrap_samples = self.n_bootstrap * (self.max_depth + 1);
        let mut bootstrap_weights: Vec<Float> = vec![0.0; n_bootstrap_samples];

        if !scene.lights.is_empty() {
            let progress = create_progress_bar(self.n_bootstrap as u64);
            progress.set_message("Generating bootstrap paths");

            thread::scope(|scope| {
                let chunk_size = clamp(self.n_bootstrap / 128, 1, 8192);
                let (tx_collector, rx_collector) = crossbeam_channel::bounded::<(usize, Float)>(chunk_size);
                let (tx_worker, rx_worker) = crossbeam_channel::bounded::<usize>(chunk_size);

                // Spawn collector thread.
                let bootstrap_weights: &mut [Float] = bootstrap_weights.as_mut_slice();
                scope.spawn(move || {
                    for (rng_index, bootstrap_weight) in rx_collector.iter() {
                        bootstrap_weights[rng_index] = bootstrap_weight;
                    }

                    // Report per thread statistics.
                    report_stats!();
                });

                // Spawn worker threads.
                for _ in 0..n_threads {
                    let rx_worker = rx_worker.clone();
                    let tx_collector = tx_collector.clone();
                    let light_distr = Arc::clone(&light_distr);
                    let light_to_index = Arc::clone(&light_to_index);
                    let progress = &progress;
                    scope.spawn(move || {
                        for i in rx_worker.iter() {
                            for depth in 0..=self.max_depth {
                                let rng_index = i * (self.max_depth + 1) + depth;
                                let mut sampler = MLTSampler::new(
                                    self.mutations_per_pixel,
                                    rng_index as u64,
                                    self.sigma,
                                    self.large_step_probability,
                                    N_SAMPLE_STREAMS,
                                );

                                let ld = Arc::clone(&light_distr);
                                let li = Arc::clone(&light_to_index);

                                let (_p_raster, spectrum) = self.l(scene, ld, li, &mut sampler, depth);
                                let bootstrap_weight = spectrum.y();

                                tx_collector.send((rng_index, bootstrap_weight)).unwrap();
                            }

                            progress.inc(1);
                        }

                        // Report per thread statistics.
                        report_stats!();
                    });
                }
                drop(rx_worker); // Drop extra since we've cloned one for each woker.
                drop(tx_collector);

                // Send work.
                for i in 0..self.n_bootstrap {
                    tx_worker.send(i).unwrap();
                }
            });

            progress.finish_with_message("Bootstrap paths complete");
        }
        let bootstrap = Arc::new(Distribution1D::new(bootstrap_weights));
        let b = bootstrap.func_int * (self.max_depth + 1) as Float;

        // Run `n_chains` Markov chains in parallel.
        let film = &self.camera.get_data().film;
        let n_total_mutations = self.mutations_per_pixel as u64 * film.get_sample_bounds().area() as u64;

        const PROGRESS_FREQUENCY: u64 = 32768;
        let progress = create_progress_bar(n_total_mutations / PROGRESS_FREQUENCY);
        progress.set_message("Rendering scene");

        if !scene.lights.is_empty() {
            thread::scope(|scope| {
                let (tx_worker, rx_worker) = crossbeam_channel::bounded::<usize>(n_threads);

                // Spawn worker threads.
                for _ in 0..n_threads {
                    let rx_worker = rx_worker.clone();
                    let light_distr = Arc::clone(&light_distr);
                    let light_to_index = Arc::clone(&light_to_index);
                    let bootstrap = Arc::clone(&bootstrap);
                    let progress = &progress;
                    scope.spawn(move || {
                        for i in rx_worker.iter() {
                            let n_chain_mutations = min(
                                (i + 1) as u64 * n_total_mutations / self.n_chains as u64,
                                n_total_mutations,
                            ) - i as u64 * n_total_mutations / self.n_chains as u64;

                            // Follow i^th Markov chain for `n_chain_mutations`.

                            // Select initial state from the set of bootstrap samples.
                            let mut rng = RNG::new(i as u64);
                            let (bootstrap_index, _, _) = bootstrap.sample_discrete(rng.uniform_float());
                            let depth = bootstrap_index % (self.max_depth + 1);

                            // Initialize local variables for selected state
                            let mut sampler = MLTSampler::new(
                                self.mutations_per_pixel,
                                bootstrap_index as u64,
                                self.sigma,
                                self.large_step_probability,
                                N_SAMPLE_STREAMS,
                            );

                            let (mut p_current, mut l_current) = self.l(
                                scene,
                                Arc::clone(&light_distr),
                                Arc::clone(&light_to_index),
                                &mut sampler,
                                depth,
                            );

                            // Run the Markov chain for `n_chain_mutations` steps.
                            for j in 0..n_chain_mutations as usize {
                                sampler.start_iteration();
                                let (p_proposed, l_proposed) = self.l(
                                    scene,
                                    Arc::clone(&light_distr),
                                    Arc::clone(&light_to_index),
                                    &mut sampler,
                                    depth,
                                );

                                // Compute acceptance probability for proposed sample.
                                let accept = min(1.0, l_proposed.y() / l_current.y());

                                // Splat both current and proposed samples to `film`.
                                if accept > 0.0 {
                                    film.add_splat(&p_proposed, &(l_proposed * accept / l_proposed.y()));
                                }
                                film.add_splat(&p_current, &(l_current * (1.0 - accept) / l_current.y()));

                                // Accept or reject the proposal.
                                if rng.uniform_float() < accept {
                                    p_current = p_proposed;
                                    l_current = l_proposed;
                                    sampler.accept();
                                } else {
                                    sampler.reject();
                                }

                                if (i * n_total_mutations as usize / self.n_chains + j) % PROGRESS_FREQUENCY as usize
                                    == 0
                                {
                                    progress.inc(1);
                                }
                            }
                        }

                        // Report per thread statistics.
                        report_stats!();
                    });
                }
                drop(rx_worker); // Drop extra since we've cloned one for each woker.

                // Send work.
                for i in 0..self.n_chains {
                    tx_worker.send(i).unwrap();
                }
            });
        }

        // Store final image computed with MLT
        progress.set_message("Writing image");
        film.write_image(b / self.mutations_per_pixel as Float);
        progress.inc(1);

        progress.finish_with_message("Render complete");
    }

    /// Preprocess the scene.
    ///
    /// * `scene` - The scene
    fn preprocess(&mut self, _scene: &Scene) {}

    /// Returns the incident radiance at the origin of a given ray.
    ///
    /// * `ray`     - The ray.
    /// * `scene`   - The scene.
    /// * `sampler` - The sampler.
    /// * `depth`   - The recursion depth.
    fn li(&self, _ray: &mut Ray, _scene: &Scene, _sampler: &mut dyn Sampler, _depth: usize) -> Spectrum {
        unimplemented!("Integrator::li() not implemented for MLTIntegrator")
    }
}

impl From<(&ParamSet, ArcCamera)> for MLTIntegrator {
    /// Create a `MLTIntegrator ` from given parameter set and camera.
    ///
    /// * `p` - A tuple containing parameter set and camera.
    fn from(p: (&ParamSet, ArcCamera)) -> Self {
        let (params, camera) = p;

        let max_depth = params.find_one_int("maxdepth", 5) as usize;
        let mut n_bootstrap = params.find_one_int("bootstrapsamples", 100000) as usize;
        let n_chains = params.find_one_int("chains", 1000) as usize;
        let mut mutations_per_pixel = params.find_one_int("mutationsperpixel", 100) as usize;
        let large_step_probability = params.find_one_float("largestepprobability", 0.3);
        let sigma = params.find_one_float("sigma", 0.01);

        if OPTIONS.quick_render {
            mutations_per_pixel = max(1, mutations_per_pixel / 16);
            n_bootstrap = max(1, n_bootstrap / 16);
        }

        MLTIntegrator::new(
            camera,
            max_depth,
            n_bootstrap,
            n_chains,
            mutations_per_pixel,
            sigma,
            large_step_probability,
        )
    }
}
