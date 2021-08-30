//! BSDF

#![allow(dead_code)]
use super::*;
use crate::rng::*;

/// Stores combinations of reflection models.
pub type BxDFType = u8;

/// Types of BSDF models.
pub const BSDF_NONE: BxDFType = 0;
pub const BSDF_REFLECTION: BxDFType = 1 << 0;
pub const BSDF_TRANSMISSION: BxDFType = 1 << 1;
pub const BSDF_DIFFUSE: BxDFType = 1 << 2;
pub const BSDF_GLOSSY: BxDFType = 1 << 3;
pub const BSDF_SPECULAR: BxDFType = 1 << 4;
pub const BSDF_ALL: BxDFType =
    BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR | BSDF_REFLECTION | BSDF_TRANSMISSION;

/// Maximum number of BxDFs that can be stored in `BSDF`.
pub const MAX_BXDFS: usize = 8;

/// BSDF interface represents a collection of BRDFs and BTDFs.
#[derive(Clone)]
pub struct BSDF {
    /// The shading normal given by per-vertex normals and/or bump mapping.
    /// It is the first axis in the orthonormal coordinate system and also
    /// used to define hemispheres for integrating incident illumincation for
    /// surface reflection.
    pub ns: Normal3f,

    /// The geometric normal defined by surface geometry.
    pub ng: Normal3f,

    /// Second axis for the orthonormal coordinate system.
    pub ss: Vector3f,

    /// Third axis for the orthonormal coordinate system.
    pub ts: Vector3f,

    /// The `BxDFs`.
    pub bxdfs: Vec<ArcBxDF>,

    /// Relative index of refraction over the surfaceboundary.
    pub eta: Float,
}

impl BSDF {
    /// Creates a new `BSDF`.
    ///
    /// * `si`  - The differential geometry at the point on a surface.
    /// * `eta` - Optional relative index of refraction over the surface
    ///           boundary. If not provided, defaults to 1.0; used for
    ///           opaque surfaces.
    pub fn new(si: &SurfaceInteraction, eta: Option<Float>) -> Self {
        let eta = eta.map_or_else(|| 1.0, |e| e);
        let ns = si.shading.n;
        let ss = si.shading.dpdu.normalize();

        Self {
            eta,
            ns,
            ng: si.hit.n,
            ss,
            ts: Vector3::from(ns).cross(&ss),
            bxdfs: Vec::with_capacity(MAX_BXDFS),
        }
    }

    /// Add a `BxDF`.
    ///
    /// * `bxdf` - The `BxDF`.
    pub fn add(&mut self, bxdf: ArcBxDF) {
        assert!(
            self.bxdfs.len() < MAX_BXDFS,
            "Cannot add BxDFs. BSDF maximum limit {} reached.",
            MAX_BXDFS
        );
        self.bxdfs.push(Arc::clone(&bxdf));
    }

    /// Returns the number of `BxDF`s that match the given type.
    ///
    /// * `bxdf_type` - The `BxdFType` to match (default to `BSDF_ALL`).
    pub fn num_components(&self, bxdf_type: BxDFType) -> usize {
        self.bxdfs
            .iter()
            .filter(|b| b.get_type() == bxdf_type)
            .count()
    }

    /// Transforms a vector from world space to local space.
    ///
    /// * `v` - The vector to transform.
    pub fn world_to_local(&self, v: &Vector3f) -> Vector3f {
        Vector3f::new(v.dot(&self.ss), v.dot(&self.ts), v.dot(&self.ns))
    }

    /// Transforms a vector from local space to world space.
    ///
    /// * `v` - The vector to transform.
    pub fn local_to_world(&self, v: &Vector3f) -> Vector3f {
        Vector3f::new(
            self.ss.x * v.x + self.ts.x * v.y + self.ns.x * v.z,
            self.ss.y * v.x + self.ts.y * v.y + self.ns.y * v.z,
            self.ss.z * v.x + self.ts.z * v.y + self.ns.z * v.z,
        )
    }

    /// Returns the BSDF evaluated for a pair of directions.
    ///
    /// * `wo_w`      - Outgoing direction in world-space.
    /// * `wi_w`      - Incident direction in world-space.
    /// * `bxdf_type` - The `BxdFType` to evaluate.
    pub fn f(&self, wo_w: &Vector3f, wi_w: &Vector3f, bxdf_type: BxDFType) -> Spectrum {
        let wi = self.world_to_local(wi_w);
        let wo = self.world_to_local(wo_w);

        if wo.z == 0.0 {
            Spectrum::new(0.0)
        } else {
            let reflect = wi_w.dot(&self.ng) * wo_w.dot(&self.ng) > 0.0;
            let mut l = Spectrum::new(0.0);
            for bxdf in self.bxdfs.iter() {
                if bxdf.matches_flags(bxdf_type)
                    && ((reflect && bxdf.get_type() & BSDF_REFLECTION > 0)
                        || (!reflect && bxdf.get_type() & BSDF_TRANSMISSION > 0))
                {
                    l += bxdf.f(&wo, &wi);
                }
            }
            l
        }
    }

    /// Returns the value of the BSDF given the outgpoing direction.
    /// direction.
    ///
    /// * `wo_w`      - Outgoing direction in world-space.
    /// * `u`         - The 2D uniform random values.
    /// * `bxdf_type` - The `BxdFType` to evaluate.
    pub fn sample_f(&self, wo_w: &Vector3f, u: &Point2f, bxdf_type: BxDFType) -> BxDFSample {
        // Choose which `BxDF` to sample.
        let matching_comps = self.num_components(bxdf_type);
        if matching_comps == 0 {
            return BxDFSample::default();
        }
        let comp = min(
            (u[0] * matching_comps as Float).floor() as usize,
            matching_comps - 1,
        );

        // Get BxDF for chosen component.
        let mut count = comp;
        let mut bxdf: Option<ArcBxDF> = None;
        for b in self.bxdfs.iter() {
            if b.matches_flags(bxdf_type) && count == 0 {
                bxdf = Some(Arc::clone(b));
                break;
            }
            count -= 1;
        }
        let bxdf = bxdf.expect("bsdf::sample_f() did not find matching bxdf");

        // Remap BxDF sample `u` to `[0,1)^2`.
        let u_remapped = Point2f::new(
            min(u[0] * (matching_comps - comp) as Float, ONE_MINUS_EPSILON),
            u[1],
        );

        // Sample chosen `BxDF`.
        let wo = self.world_to_local(wo_w);
        if wo.z == 0.0 {
            return BxDFSample::default();
        }

        let sampled_type = bxdf.get_type();
        let sample = bxdf.sample_f(&wo, &u_remapped);
        let mut pdf = sample.pdf;
        if pdf == 0.0 {
            return BxDFSample::from(sampled_type);
        }
        let wi_world = self.local_to_world(&sample.wi);

        // Compute overall PDF with all matching BxDFs.
        if bxdf.get_type() & BSDF_SPECULAR == 0 && matching_comps > 1 {
            for b in self.bxdfs.iter() {
                if Arc::ptr_eq(&b, &bxdf) && b.matches_flags(bxdf_type) {
                    pdf += b.pdf(&wo, &sample.wi);
                }
            }
        }
        if matching_comps > 1 {
            pdf /= matching_comps as Float;
        }

        // Compute value of BSDF for sampled direction.
        let mut f = Spectrum::new(0.0);
        if bxdf.get_type() & BSDF_SPECULAR == 0 {
            let reflect = wi_world.dot(&self.ng) * wo_w.dot(&self.ng) > 0.0;
            for bxdf in self.bxdfs.iter() {
                if bxdf.matches_flags(bxdf_type)
                    && ((reflect && bxdf.get_type() & BSDF_REFLECTION > 0)
                        || (!reflect && bxdf.get_type() & BSDF_TRANSMISSION > 0))
                {
                    f += bxdf.f(&wo, &sample.wi);
                }
            }
        }
        BxDFSample::new(f, pdf, wi_world, sampled_type)
    }

    /// Computes the hemispherical-directional reflectance function ρ.
    ///
    /// * `wo_w`      - Outgoing direction in world-space.
    /// * `u`         - Samples used by Monte Carlo algorithm.
    /// * `bxdf_type` - The `BxdFType` to evaluate.
    pub fn rho_hd(&self, wo_w: &Vector3f, u: &[Point2f], bxdf_type: BxDFType) -> Spectrum {
        let wo = self.world_to_local(wo_w);

        let mut l = Spectrum::new(0.0);
        for bxdf in self.bxdfs.iter() {
            if bxdf.matches_flags(bxdf_type) {
                l += bxdf.rho_hd(&wo, u);
            }
        }
        l
    }

    /// Computes the hemispherical-hemispherical-directional reflectance function ρ.
    ///
    /// * `u1`        - Samples used b Monte Carlo algorithm.
    /// * `u2`        - Samples used b Monte Carlo algorithm.
    /// * `bxdf_type` - The `BxdFType` to evaluate.
    pub fn rho_hh(&self, u1: &[Point2f], u2: &[Point2f], bxdf_type: BxDFType) -> Spectrum {
        let mut l = Spectrum::new(0.0);
        for bxdf in self.bxdfs.iter() {
            if bxdf.matches_flags(bxdf_type) {
                l += bxdf.rho_hh(u1, u2);
            }
        }
        l
    }

    /// Evaluates the PDF for the sampling method. Default is based on the
    /// cosine-weighted sampling in `BxDF::sample_f()` default implementation.
    ///
    /// * `wo_w`      - Outgoing direction in world-space.
    /// * `wi_w`      - Incident direction in world-space.
    /// * `bxdf_type` - The `BxdFType` to evaluate.
    pub fn pdf(&self, wo_w: &Vector3f, wi_w: &Vector3f, bxdf_type: BxDFType) -> Float {
        if self.bxdfs.len() == 0 {
            return 0.0;
        }

        let wo = self.world_to_local(wo_w);
        let wi = self.world_to_local(wi_w);

        if wo.z == 0.0 {
            return 0.0;
        }

        let mut matching_comps = 0;
        let mut pdf = 0.0;
        for bxdf in self.bxdfs.iter() {
            if bxdf.matches_flags(bxdf_type) {
                matching_comps += 1;
                pdf += bxdf.pdf(&wo, &wi);
            }
        }
        if matching_comps > 0 {
            pdf / matching_comps as Float
        } else {
            0.0
        }
    }
}

/// Atomic reference counted `BSDF`.
pub type ArcBSDF = Arc<BSDF>;
