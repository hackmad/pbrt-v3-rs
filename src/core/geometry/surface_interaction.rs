//! Surface Interactions

#![allow(dead_code)]
use crate::core::bssrdf::*;
use crate::core::geometry::*;
use crate::core::material::*;
use crate::core::pbrt::*;
use crate::core::primitive::*;
use crate::core::reflection::*;
use crate::core::spectrum::*;

/// SurfaceInteraction represents geometry of a particular point on a surface.
#[derive(Clone)]
pub struct SurfaceInteraction<'a> {
    /// The common interaction data.
    pub hit: Hit,

    /// The uv coordinates from surface parametrization.
    pub uv: Point2f,

    /// Parametric partial derivative of the point ∂p/∂u.
    pub dpdu: Vector3f,

    /// Parametric partial derivative of the point ∂p/∂v.
    pub dpdv: Vector3f,

    /// Differential change ∂n/∂u in surface normal as we move along u.
    pub dndu: Normal3f,

    /// Differential change ∂n/∂v in surface normal as we move along v.
    pub dndv: Normal3f,

    /// Differential change ∂u/∂x in parameteric coordinate u as we move along x.
    pub dudx: Float,

    /// Differential change ∂u/∂y in parameteric coordinate u as we move along y.
    pub dudy: Float,

    /// Differential change ∂v/∂x in parameteric coordinate v as we move along x.
    pub dvdx: Float,

    /// Differential change ∂v/∂y in parameteric coordinate v as we move along y.
    pub dvdy: Float,

    /// Partial derivative of the point ∂p/∂x in world space.
    pub dpdx: Vector3f,

    /// Partial derivative of the point ∂p/∂y in world space.
    pub dpdy: Vector3f,

    /// Shading geometry used for perturbed values.
    pub shading: Shading,

    /// The shape.
    pub shape: Option<ArcShape>,

    /// The BSDF.
    pub bsdf: Option<ArcBSDF>,

    /// The BSSRDF.
    pub bssrdf: Option<ArcBSSRDF>,

    /// The primitive.
    pub primitive: Option<&'a dyn Primitive>,
}

impl<'a> SurfaceInteraction<'a> {
    /// Create a new surface interaction.
    ///
    /// * `p`                - Point of interaction.
    /// * `p_error`          - Floating point error for ray intersection points.
    /// * `uv`               - The uv coordinates from surface parametrization.
    /// * `wo`               - The negative ray direction (outgoing direction used
    ///                        when computing lighting at points).
    /// * `dpdu`             - Parametric partial derivative of the point ∂p/∂u.
    /// * `dpdv`             - Parametric partial derivative of the point ∂p/∂v.
    /// * `dndu`             - Differential change ∂n/∂v in surface normal as we move along u.
    /// * `dndv`             - Differential change ∂n/∂v in surface normal as we move along v.
    /// * `time`             - Time when interaction occurred.
    /// * `shape`            - The shape.
    pub fn new(
        p: Point3f,
        p_error: Vector3f,
        uv: Point2f,
        wo: Vector3f,
        dpdu: Vector3f,
        dpdv: Vector3f,
        dndu: Normal3f,
        dndv: Normal3f,
        time: Float,
        shape: Option<ArcShape>,
    ) -> Self {
        // Calculate normal n from the partial derivatives.
        let mut n = Normal3f::from(dpdu.cross(&dpdv).normalize());

        // Adjust normal based on orientation and handedness
        if let Some(s) = shape.clone() {
            if s.get_data().reverse_orientation ^ s.get_data().transform_swaps_handedness {
                n *= -1.0;
            }
        }

        Self {
            hit: Hit::new(p, time, p_error, wo, n, None),
            uv,
            dpdu,
            dpdv,
            dndu,
            dndv,
            dudx: 0.0,
            dudy: 0.0,
            dvdx: 0.0,
            dvdy: 0.0,
            dpdx: Vector3f::default(),
            dpdy: Vector3f::default(),
            shape: shape.clone(),
            shading: Shading::new(n, dpdu, dpdv, dndu, dndv),
            bsdf: None,
            bssrdf: None,
            primitive: None,
        }
    }

    /// Returns updated shading geometry.
    ///
    /// * `dpdu` - Parametric partial derivative of the point ∂p/∂u.
    /// * `dpdv` - Parametric partial derivative of the point ∂p/∂v.
    /// * `dndu` - Differential change ∂n/∂v in surface normal as we move along u.
    /// * `dndv` - Differential change ∂n/∂v in surface normal as we move along v.
    pub fn set_shading_geometry(
        &mut self,
        dpdu: Vector3f,
        dpdv: Vector3f,
        dndu: Normal3f,
        dndv: Normal3f,
        orientation_is_authoritative: bool,
    ) {
        // Compute normal.
        let mut hit_n = self.hit.n;
        let mut shading_n = Normal3::from(dpdu.cross(&dpdv)).normalize();

        if let Some(s) = self.shape.clone() {
            if s.get_data().reverse_orientation ^ s.get_data().transform_swaps_handedness {
                shading_n = -self.shading.n;
                if orientation_is_authoritative {
                    hit_n = hit_n.face_forward(&shading_n.into());
                } else {
                    shading_n = shading_n.face_forward(&hit_n.into());
                }
            }
        }

        // Initialize shading partial derivative values.
        self.hit.n = hit_n;
        self.shading = Shading::new(shading_n, dpdu, dpdv, dndu, dndv);
    }

    /// Initializes representations of the light-scattering properties of the
    /// material at the intersection point on the primtive's surface.
    ///
    /// * `ray`                  - The ray.
    /// * `mode`                 - Transport mode.
    /// * `allow_multiple_lobes` - Indicates whether the material should use
    ///                            BxDFs that aggregate multiple types of
    ///                            scattering into a single BxDF when such BxDFs
    ///                            are available.
    pub fn compute_scattering_functions(
        &mut self,
        ray: &Ray,
        allow_multiple_lobes: bool,
        mode: TransportMode,
    ) {
        self.compute_differentials(ray);
        if let Some(primitive) = self.primitive {
            primitive.compute_scattering_functions(self, mode, allow_multiple_lobes);
        }
    }

    /// Use offset rays to estimate the partial derivatives mapping p(x, y) from
    /// image position to world space position and the partial derivatives of the
    /// mappings u(x, y) and v(x, y) from (x, y) to (u, v) parametric coordinates,
    /// giving theworld space positions ∂p/∂x and ∂p/∂y.
    ///
    /// * `ray` - The ray.
    pub fn compute_differentials(&mut self, ray: &Ray) {
        if let Some(rd) = ray.differentials {
            // Estimate screen space change in `pt{}` and `(u,v)`.
            let n = self.hit.n;
            let p = self.hit.p;

            // Compute auxiliary intersection points with plane.
            let d = self.hit.n.dot(&Vector3f::new(p.x, p.y, p.z));
            let tx = -(n.dot(&Vector3f::from(rd.rx_origin)) - d) / n.dot(&rd.rx_direction);
            if tx.is_infinite() || tx.is_nan() {
                self.dudx = 0.0;
                self.dvdx = 0.0;
                self.dudy = 0.0;
                self.dvdy = 0.0;
                self.dpdx = Vector3f::new(0.0, 0.0, 0.0);
                self.dpdy = Vector3f::new(0.0, 0.0, 0.0);
                return;
            }
            let px = rd.rx_origin + tx * rd.rx_direction;

            let ty = -(n.dot(&Vector3f::from(rd.ry_origin)) - d) / n.dot(&rd.ry_direction);
            if ty.is_infinite() || ty.is_nan() {
                self.dudx = 0.0;
                self.dvdx = 0.0;
                self.dudy = 0.0;
                self.dvdy = 0.0;
                self.dpdx = Vector3f::new(0.0, 0.0, 0.0);
                self.dpdy = Vector3f::new(0.0, 0.0, 0.0);
                return;
            }
            let py = rd.ry_origin + ty * rd.ry_direction;

            self.dpdx = px - p;
            self.dpdy = py - p;

            // Compute `(u,v)` offsets at auxiliary points.

            // Choose two dimensions to use for ray offset computation.
            let dim = if abs(n.x) > abs(n.y) && abs(n.x) > abs(n.z) {
                [1, 2]
            } else if abs(n.y) > abs(n.z) {
                [0, 2]
            } else {
                [0, 1]
            };

            // Initialize `A`, `Bx`, and `By` matrices for offset computation.
            let a = [
                [self.dpdu[dim[0]], self.dpdv[dim[0]]],
                [self.dpdu[dim[1]], self.dpdv[dim[1]]],
            ];
            let bx = [px[dim[0]] - p[dim[0]], px[dim[1]] - p[dim[1]]];
            let by = [py[dim[0]] - p[dim[0]], py[dim[1]] - p[dim[1]]];
            if let Some((dudx, dvdx)) = solve_linear_system_2x2(&a, &bx) {
                self.dudx = dudx;
                self.dvdx = dvdx;
            } else {
                self.dudx = 0.0;
                self.dvdx = 0.0;
            }
            if let Some((dudy, dvdy)) = solve_linear_system_2x2(&a, &by) {
                self.dudy = dudy;
                self.dvdy = dvdy;
            } else {
                self.dudy = 0.0;
                self.dvdy = 0.0;
            }
        } else {
            self.dudx = 0.0;
            self.dvdx = 0.0;
            self.dudy = 0.0;
            self.dvdy = 0.0;
            self.dpdx = Vector3f::new(0.0, 0.0, 0.0);
            self.dpdy = Vector3f::new(0.0, 0.0, 0.0);
        }
    }

    /// Returns the emitted radiance at a surface point intersected by a ray
    /// for an area light.
    ///
    /// * `w` - The outgoing direction.
    pub fn le(&self, w: &Vector3f) -> Spectrum {
        if let Some(area_light) = self.primitive.map(|p| p.get_area_light()).flatten() {
            area_light.l(self.get_hit(), &w)
        } else {
            Spectrum::new(0.0)
        }
    }
}

impl<'a> Interaction for SurfaceInteraction<'a> {
    /// Returns the interaction hit data.
    fn get_hit(&self) -> &Hit {
        &self.hit
    }

    /// Returns the surface interaction.
    ///
    /// NOTE: This is a hack because I don't want to write a function
    /// for retrieving every field in `SurfaceInteraction`. If there is a clean
    /// way of retrieving a struct that implements an interface I will get rid
    /// of this.
    fn get_surface_interaction(&self) -> Option<&SurfaceInteraction> {
        Some(self)
    }

    /// Returns the medium interaction or None.
    ///
    /// NOTE: This is a hack because I don't want to write a function
    /// for retrieving every field in `MediumInteraction`. If there is a clean
    /// way of retrieving a struct that implements an interface I will get rid
    /// of this.
    fn get_medium_interaction(&self) -> Option<&MediumInteraction> {
        None
    }
}

/// Shading geometry used for perturbed values for bump mapping.
#[derive(Clone)]
pub struct Shading {
    /// Surface normal.
    pub n: Normal3f,

    /// Parametric partial derivative of the point ∂p/∂u.
    pub dpdu: Vector3f,

    /// Parametric partial derivative of the point ∂p/∂v.
    pub dpdv: Vector3f,

    /// Differential change ∂n/∂v in surface normal as we move along u.
    pub dndu: Normal3f,

    /// Differential change ∂n/∂v in surface normal as we move along v.
    pub dndv: Normal3f,
}

impl Shading {
    /// Create a new shading struct.
    /// * `n`    - Surface normal.
    /// * `dpdu` - Parametric partial derivative of the point ∂p/∂u.
    /// * `dpdv` - Parametric partial derivative of the point ∂p/∂v.
    /// * `dndu` - Differential change ∂n/∂v in surface normal as we move along u.
    /// * `dndv` - Differential change ∂n/∂v in surface normal as we move along v.
    pub fn new(
        n: Normal3f,
        dpdu: Vector3f,
        dpdv: Vector3f,
        dndu: Normal3f,
        dndv: Normal3f,
    ) -> Self {
        Self {
            n,
            dpdu,
            dpdv,
            dndu,
            dndv,
        }
    }
}
