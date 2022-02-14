//! Bidirectional scattering surface reflectance distribution function.

use crate::geometry::*;
use crate::interaction::*;
use crate::material::*;
use crate::pbrt::*;
use crate::reflection::*;
use crate::scene::Scene;
use crate::spectrum::*;
use bumpalo::Bump;
use std::sync::Arc;

/// Implements a simple BSSRDF implementation that can support general `Shapes`.
pub struct SeparableBSSRDF {
    /// Current outgoing surface interaction hit point.
    pub po_hit: Hit,

    /// Current outgoing surface interaction shading information.
    pub po_shading: Shading,

    /// Index of refraction of the scattering medium.
    pub eta: Float,

    /// Local coordinate frame: normal at hit point.
    pub ns: Normal3f,

    /// Local coordinate frame: shading parametric partial derivative of the
    /// point ∂p/∂u (normalized).
    pub ss: Vector3f,

    /// Local coordinate frame: cross product of `ns` x `ss`.
    pub ts: Vector3f,

    /// The underlying material.
    pub material: ArcMaterial,

    /// Light transport mode.
    pub mode: TransportMode,
}

impl SeparableBSSRDF {
    /// Allocate a new instance of `SeparableBSSRDF`.
    ///
    /// * `po`       - Current outgoing surface interaction.
    /// * `eta`      - Index of refraction of the scattering medium.
    /// * `material` - The material.
    /// * `mode`     - Light transport mode.
    pub fn alloc<'arena>(
        arena: &'arena Bump,
        po: &SurfaceInteraction,
        eta: Float,
        material: ArcMaterial,
        mode: TransportMode,
    ) -> &'arena mut Self {
        let ns = po.shading.n;
        let ss = po.shading.dpdu.normalize();
        let ts = Vector3f::from(&ns.cross(&ss));

        arena.alloc(Self {
            po_hit: po.hit.clone(),
            po_shading: po.shading.clone(),
            eta,
            ns,
            ss,
            ts,
            material,
            mode,
        })
    }

    /// Clone into a newly allocated a new instance of `TabulatedBSSRDF`.
    ///
    /// * `arena` - The arena for memory allocations.
    #[allow(clippy::mut_from_ref)]
    pub fn clone_alloc<'arena>(&self, arena: &'arena Bump) -> &'arena mut Self {
        arena.alloc(Self {
            po_hit: self.po_hit.clone(),
            po_shading: self.po_shading.clone(),
            eta: self.eta,
            ns: self.ns.clone(),
            ss: self.ss.clone(),
            ts: self.ts.clone(),
            material: Arc::clone(&self.material),
            mode: self.mode,
        })
    }

    /// Evaluates the eight-dimensional distribution function S(), which quantifies / the ratio of differential radiance at point `po` in direction `ωo` to the
    /// incident differential flux at `pi` from direction `ωi`.
    ///
    /// * `pi` - Interaction point for incident differential flux.
    /// * `wi` - Direction for incident different flux.
    /// * `sr` - Evaluates the radial profile function based on distance between points.
    pub fn s<Sr>(&self, pi: &SurfaceInteraction, wi: &Vector3f, sr: Sr) -> Spectrum
    where
        Sr: Fn(Float) -> Spectrum,
    {
        let ft = 1.0 - fr_dielectric(self.po_hit.wo.dot(&self.po_shading.n), 1.0, self.eta);
        ft * self.sp(pi, sr) * self.sw(wi)
    }

    /// Evaluates the spatial term for the distribution function.
    ///
    /// * `pi` - Interaction point for incident differential flux.
    /// * `sr` - Evaluates the radial profile function based on distance between points.
    pub fn sp<Sr>(&self, pi: &SurfaceInteraction, sr: Sr) -> Spectrum
    where
        Sr: Fn(Float) -> Spectrum,
    {
        sr(self.po_hit.p.distance(pi.hit.p))
    }

    /// Evaluates the directional term for the distribution function.
    ///
    /// * `w` - Direction for incident different flux.
    pub fn sw(&self, w: &Vector3f) -> Spectrum {
        let c = 1.0 - 2.0 * fresnel_moment_1(1.0 / self.eta);
        Spectrum::new((1.0 - fr_dielectric(cos_theta(w), 1.0, self.eta)) / (c * PI))
    }

    /// Returns the value of the BSSRDF, the surface position where a ray
    /// re-emerges following internal scattering and probability density function.
    ///
    /// * `scene`     - The scene.
    /// * `u1`        - Sample values for Monte Carlo.
    /// * `u2`        - Sample values for Monte Carlo.
    /// * `si`        - The surface position where a ray re-emerges following internal
    ///                 scattering.
    /// * `sample_sr` - Samples radius values proportional to the radial profile
    ///                 function.
    /// * `pdf_sp`    - Evaluate the combined PDF that takes all of the sampling
    ///                 strategies `sample_sp()` into account.
    /// * `sr`        - Evaluates the radial profile function based on distance
    ///                 between points.
    pub fn sample_s<'scene, SampleSr, PdfSp, Sr>(
        &self,
        scene: &'scene Scene,
        u1: Float,
        u2: &Point2f,
        si: &'scene mut SurfaceInteraction,
        sample_sr: SampleSr,
        pdf_sp: PdfSp,
        sr: Sr,
    ) -> (Spectrum, Float)
    where
        SampleSr: Fn(usize, Float) -> Float,
        PdfSp: Fn(&SurfaceInteraction) -> Float,
        Sr: Fn(Float) -> Spectrum,
    {
        let (sp, pdf) = self.sample_sp(scene, u1, u2, si, sample_sr, pdf_sp, sr);

        if !sp.is_black() {
            // Initialize material model at sampled surface interaction.
            //let bsdf = Some(BSDF::alloc(arena, si, None));
            //si.bsdf = bsdf;
            si.hit.wo = Vector3f::from(si.shading.n);
        }

        (sp, pdf)
    }

    /// Use a different sampling technique per wavelength to deal with spectral
    /// variation, and each technique is additionally replicated three times with
    /// different projection axes given by the basis vectors of a local frame,
    /// resulting in a total of 3 * Spectrum::nSamples sampling techniques. This
    /// ensures that every point `S` where takes on non-negligible values is
    /// intersected with a reasonable probability.
    ///
    /// * `scene`     - The scene.
    /// * `u1`        - Sample values for Monte Carlo.
    /// * `u2`        - Sample values for Monte Carlo.
    /// * `si`        - Surface interaction representing sample returned.
    /// * `sample_sr` - Samples radius values proportional to the radial profile
    ///                 function.
    /// * `pdf_sp`    - Evaluate the combined PDF that takes all of the sampling
    ///                 strategies `sample_sp()` into account.
    /// * `sr`        - Evaluates the radial profile function based on distance
    ///                 between points.
    pub fn sample_sp<'scene, SampleSr, PdfSp, Sr>(
        &self,
        scene: &'scene Scene,
        u1: Float,
        u2: &Point2f,
        si: &'scene mut SurfaceInteraction<'_, '_>,
        sample_sr: SampleSr,
        pdf_sp: PdfSp,
        sr: Sr,
    ) -> (Spectrum, Float)
    where
        SampleSr: Fn(usize, Float) -> Float,
        PdfSp: Fn(&SurfaceInteraction) -> Float,
        Sr: Fn(Float) -> Spectrum,
    {
        // Choose projection axis for BSSRDF sampling.
        let (vx, vy, vz, mut u1) = if u1 < 0.5 {
            (self.ss, self.ts, Vector3f::from(self.ns), u1 * 2.0)
        } else if u1 < 0.75 {
            // Prepare for sampling rays with respect to `ss`.
            (self.ts, Vector3f::from(self.ns), self.ss, (u1 - 0.5) * 4.0)
        } else {
            // Prepare for sampling rays with respect to `ts`.
            (Vector3f::from(self.ns), self.ss, self.ts, (u1 - 0.75) * 4.0)
        };

        // Choose spectral channel for BSSRDF sampling.
        let ch = clamp(
            (u1 * SPECTRUM_SAMPLES as Float) as usize,
            0,
            SPECTRUM_SAMPLES - 1,
        );
        u1 = u1 * (SPECTRUM_SAMPLES - ch) as Float;

        // Sample BSSRDF profile in polar coordinates.
        let r = sample_sr(ch, u2[0]);
        if r < 0.0 {
            return (Spectrum::ZERO, 0.0);
        }
        let phi = TWO_PI * u2[1];

        // Compute BSSRDF profile bounds and intersection height.
        let r_max = sample_sr(ch, 0.999);
        if r >= r_max {
            return (Spectrum::ZERO, 0.0);
        }
        let l = 2.0 * (r_max * r_max - r * r).sqrt();

        // Compute BSSRDF sampling ray segment
        let base_p = self.po_hit.p + r * (vx * cos(phi) + vy * sin(phi)) - l * vz * 0.5;
        let base_time = self.po_hit.time;
        let mut base = Hit::new(
            base_p,
            base_time,
            Vector3f::ZERO,
            Vector3f::ZERO,
            Normal3f::ZERO,
            None,
        );
        let p_target = base.p + l * vz;

        // Intersect BSSRDF sampling ray against the scene geometry.

        // Declare `IntersectionChain` and linked list.
        let mut chain = Vec::new();

        // Accumulate chain of intersections along ray.
        let mut n_found = 0;
        loop {
            let mut r = base.spawn_ray_to_point(&p_target);
            if r.d == Vector3f::ZERO {
                break;
            }
            let isect = scene.intersect(&mut r);
            if let Some(it) = isect {
                base = it.hit.clone();
                // Append admissible intersection to `IntersectionChain`.
                if let Some(material) = it.primitive.map(|p| p.get_material()).flatten() {
                    if Arc::ptr_eq(&material, &self.material) {
                        chain.push(it);
                        n_found += 1;
                    }
                }
            } else {
                break;
            }
        }

        // Randomly choose one of several intersections during BSSRDF sampling.
        if n_found == 0 {
            return (Spectrum::ZERO, 0.0);
        }

        let mut selected = clamp((u1 * n_found as Float) as Int, 0, n_found as Int - 1);
        let mut idx = 0;
        while selected > 0 {
            selected -= 1;
            idx += 1;
        }
        *si = chain[idx];
        let pi = si;

        // Compute sample PDF and return the spatial BSSRDF term $\Sp$.
        let pdf = pdf_sp(pi) / n_found as Float;
        let term = self.sp(pi, sr);

        (term, pdf)
    }

    /// Evaluate the combined PDF that takes all of the sampling strategies
    /// `sample_sp()` into account.
    ///
    /// * `si` - Surface interaction.
    pub fn pdf_sp(&self, pi: &SurfaceInteraction) -> Float {
        todo!()
    }
}

/// Evaluate first moment of the Fresnel reflectance function.
///
/// * `eta` - Index of refraction of the scattering medium.
pub fn fresnel_moment_1(eta: Float) -> Float {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let eta5 = eta4 * eta;
    if eta < 1.0 {
        0.45966 - 1.73965 * eta + 3.37668 * eta2 - 3.904945 * eta3 + 2.49277 * eta4 - 0.68441 * eta5
    } else {
        -4.61686 + 11.1136 * eta - 10.4646 * eta2 + 5.11455 * eta3 - 1.27198 * eta4 + 0.12746 * eta5
    }
}

/// Evaluate second moment of the Fresnel reflectance function.
///
/// * `eta` - Index of refraction of the scattering medium.
pub fn fresnel_moment_2(eta: Float) -> Float {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let eta5 = eta4 * eta;
    if eta < 1.0 {
        0.27614 - 0.87350 * eta + 1.12077 * eta2 - 0.65095 * eta3 + 0.07883 * eta4 + 0.04860 * eta5
    } else {
        let r_eta = 1.0 / eta;
        let r_eta2 = r_eta * r_eta;
        let r_eta3 = r_eta2 * r_eta;
        -547.033 + 45.3087 * r_eta3 - 218.725 * r_eta2 + 458.843 * r_eta + 404.557 * eta
            - 189.519 * eta2
            + 54.9327 * eta3
            - 9.00603 * eta4
            + 0.63942 * eta5
    }
}
