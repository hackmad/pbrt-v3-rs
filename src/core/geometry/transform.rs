//! Transformations

#![allow(dead_code)]
use super::{
    abs, gamma, matrix4x4, normal3, point3, ray, ray_differential, ray_with_differentials, shading,
    surface_interaction, vector3, Bounds3, Bounds3f, Dot, FaceForward, Float, Matrix4x4, Normal3f,
    Point3f, Ray, SurfaceInteraction, Union, Vector3, Vector3f, IDENTITY_MATRIX,
};
use std::cmp::{Eq, Ord, Ordering, PartialOrd};
use std::ops::Mul;
use std::sync::Arc;

/// A transformation for mapping from points to points and vectors to vectors.
#[derive(Copy, Clone, Debug, Default)]
pub struct Transform {
    /// The transformation matrix.
    pub m: Matrix4x4,

    /// The inverse transformation matrix.
    pub m_inv: Matrix4x4,
}

/// Atomic reference counted `Transform`.
pub type ArcTransform = Arc<Transform>;

/// Create a transformation from a 2-dimensional array representing a 4x4 matrix.
///
/// * `mat` - A matrix representing a transformation.
#[rustfmt::skip]
pub fn transform(mat: [[Float; 4]; 4]) -> Transform {
    let m = matrix4x4(
        mat[0][0], mat[0][1], mat[0][2], mat[0][3], 
        mat[1][0], mat[1][1], mat[1][2], mat[1][3],
        mat[2][0], mat[2][1], mat[2][2], mat[2][3], 
        mat[3][0], mat[3][1], mat[3][2], mat[3][3],
    );

    Transform {
        m,
        m_inv: m.inverse(),
    }
}

/// Create a transformation representing a translation.
///
/// * `delta` -  Translation.
#[rustfmt::skip]
pub fn translate(delta: &Vector3f) -> Transform {
    Transform {
        m: matrix4x4(
            1.0, 0.0, 0.0, delta.x,
            0.0, 1.0, 0.0, delta.y, 
            0.0, 0.0, 1.0, delta.z, 
            0.0, 0.0, 0.0, 1.0,
        ),
        m_inv: matrix4x4(
            1.0, 0.0, 0.0, -delta.x, 
            0.0, 1.0, 0.0, -delta.y, 
            0.0, 0.0, 1.0, -delta.z,
            0.0, 0.0, 0.0,  1.0,
        ),
    }
}

/// Create a transformation representing a scale.
///
/// * `x` -  Scaling factor in x-axis.
/// * `y` -  Scaling factor in y-axis.
/// * `z` -  Scaling factor in z-axis.
#[rustfmt::skip]
pub fn scale(x: Float, y: Float, z: Float) -> Transform {
    Transform {
        m: matrix4x4(
            x,   0.0, 0.0, 0.0, 
            0.0, y,   0.0, 0.0, 
            0.0, 0.0, z,   0.0, 
            0.0, 0.0, 0.0, 1.0,
        ),
        m_inv: matrix4x4(
            1.0 / x, 0.0,     0.0,     0.0,
            0.0,     1.0 / y, 0.0,     0.0,
            0.0,     0.0,     1.0 / z, 0.0,
            0.0,     0.0,     0.0,     1.0,
        ),
    }
}

/// Create a transformation representing rotation about the x-axis.
///
/// * `theta` -  Angle in degrees.
#[rustfmt::skip]
pub fn rotate_x(theta: Float) -> Transform {
    let r = theta.to_radians();
    let sin_theta = r.sin();
    let cos_theta = r.cos();
    let m = matrix4x4(
        1.0, 0.0,        0.0,       0.0, 
        0.0, cos_theta, -sin_theta, 0.0,
        0.0, sin_theta,  cos_theta, 0.0, 
        0.0, 0.0,        0.0,       1.0,
    );
    Transform { m, m_inv: m.transpose() }
}

/// Create a transformation representing rotation about the y-axis.
///
/// * `theta` -  Angle in degrees.
#[rustfmt::skip]
pub fn rotate_y(theta: Float) -> Transform {
    let r = theta.to_radians();
    let sin_theta = r.sin();
    let cos_theta = r.cos();
    let m = matrix4x4(
         cos_theta, 0.0, sin_theta, 0.0, 
         0.0,       1.0, 0.0,       0.0, 
        -sin_theta, 0.0, cos_theta, 0.0,
         0.0,       0.0, 0.0,       1.0,
    );
    Transform { m, m_inv: m.transpose() }
}

/// Create a transformation representing rotation about the z-axis.
///
/// * `theta` -  Angle in degrees.
#[rustfmt::skip]
pub fn rotate_z(theta: Float) -> Transform {
    let r = theta.to_radians();
    let sin_theta = r.sin();
    let cos_theta = r.cos();
    let m = matrix4x4(
        cos_theta, -sin_theta, 0.0, 0.0, 
        sin_theta,  cos_theta, 0.0, 0.0,
        0.0,        0.0,       1.0, 0.0, 
        0.0,        0.0,       0.0, 1.0,
    );
    Transform { m, m_inv: m.transpose() }
}

/// Create a transformation representing rotation about a vector.
///
/// * `theta` - Angle in degrees.
/// * `a`     - Vector.
pub fn rotate_axis(theta: Float, a: &Vector3f) -> Transform {
    let r = theta.to_radians();
    let sin_theta = r.sin();
    let cos_theta = r.cos();
    let mut m = Matrix4x4::default();

    // Compute rotation of first basis vector
    m.m[0][0] = a.x * a.x + (1.0 - a.x * a.x) * cos_theta;
    m.m[0][1] = a.x * a.y * (1.0 - cos_theta) - a.z * sin_theta;
    m.m[0][2] = a.x * a.z * (1.0 - cos_theta) + a.y * sin_theta;
    m.m[0][3] = 0.0;

    // Compute rotations of second and third basis vectors
    m.m[1][0] = a.x * a.y * (1.0 - cos_theta) + a.z * sin_theta;
    m.m[1][1] = a.y * a.y + (1.0 - a.y * a.y) * cos_theta;
    m.m[1][2] = a.y * a.z * (1.0 - cos_theta) - a.x * sin_theta;
    m.m[1][3] = 0.0;

    m.m[2][0] = a.x * a.z * (1.0 - cos_theta) - a.y * sin_theta;
    m.m[2][1] = a.y * a.z * (1.0 - cos_theta) + a.x * sin_theta;
    m.m[2][2] = a.z * a.z + (1.0 - a.z * a.z) * cos_theta;
    m.m[2][3] = 0.0;

    Transform {
        m,
        m_inv: m.transpose(),
    }
}

/// Generate a transformation to point a camera to a desired location.
///
/// * `pos`  - Position of camera.
/// * `look` - Position to point towards.
/// * `up`   - Used to orient the camera's viewing direction implied by `pos`
///            and `look`.
#[rustfmt::skip]
pub fn look_at(pos: &Point3f, look: &Point3f, up: &Vector3f) -> Transform {
    let dir = (*look - *pos).normalize();
    let right = up.normalize().cross(&dir).normalize();
    let new_up = dir.cross(&right);

    let camera_to_world = matrix4x4(
        right.x, new_up.x, dir.x, pos.x,
        right.y, new_up.y, dir.y, pos.y,
        right.z, new_up.z, dir.z, pos.z,
        0.0,     0.0,      0.0,   1.0,
    );

    Transform {
        m: camera_to_world.inverse(),
        m_inv: camera_to_world,
    }
}

// Returns true if x is in (0.999, 1.001) range as being close to 1.0.
//
// * `x` - The value to check
fn not_one(x: Float) -> bool {
    x < 0.999 || x > 1.001
}

impl Transform {
    // Returns the inverse transformation.
    pub fn inverse(&self) -> Transform {
        Transform {
            m: self.m_inv,
            m_inv: self.m,
        }
    }

    // Returns a transformation with the matrices transposed.
    pub fn transpose(&self) -> Transform {
        Transform {
            m: self.m.transpose(),
            m_inv: self.m_inv.transpose(),
        }
    }

    // Returns true if matrix is identity matrix
    pub fn is_identity(&self) -> bool {
        self.m == IDENTITY_MATRIX
    }

    /// Checks if transformation has a scaling term by transforming the
    /// coordinate axes and checking that their transformed length is close
    /// to 1.0.
    pub fn has_scale(&self) -> bool {
        let la2 = self
            .transform_vector(&vector3(1.0, 0.0, 0.0))
            .length_squared();
        let lb2 = self
            .transform_vector(&vector3(0.0, 1.0, 0.0))
            .length_squared();
        let lc2 = self
            .transform_vector(&vector3(0.0, 0.0, 1.0))
            .length_squared();

        not_one(la2) || not_one(lb2) || not_one(lc2)
    }

    /// Applies transformation to a given point.
    ///
    /// * `p` - The point.
    pub fn transform_point(&self, p: &Point3f) -> Point3f {
        let m = &self.m;
        let xp = m[0][0] * p.x + m[0][1] * p.y + m[0][2] * p.z + m[0][3];
        let yp = m[1][0] * p.x + m[1][1] * p.y + m[1][2] * p.z + m[1][3];
        let zp = m[2][0] * p.x + m[2][1] * p.y + m[2][2] * p.z + m[2][3];
        let wp = m[3][0] * p.x + m[3][1] * p.y + m[3][2] * p.z + m[3][3];

        debug_assert!(wp != 0.0, "Transformation<Point3f>: wp is not zero");

        if wp == 1.0 {
            point3(xp, yp, zp)
        } else {
            point3(xp, yp, zp) / wp
        }
    }

    /// Returns the transformed point and absolute error due to applying the
    /// transformation to a point.
    ///
    /// * `p` - The point.
    pub fn transform_point_with_error(&self, p: &Point3f) -> (Point3f, Vector3f) {
        // Compute absolute error for transformed point
        let m = &self.m;

        let x_abs_sum = abs(m[0][0] * p.x) + abs(m[0][1] * p.y) + abs(m[0][2] * p.z) + abs(m[0][3]);
        let y_abs_sum = abs(m[1][0] * p.x) + abs(m[1][1] * p.y) + abs(m[1][2] * p.z) + abs(m[1][3]);
        let z_abs_sum = abs(m[2][0] * p.x) + abs(m[2][1] * p.y) + abs(m[2][2] * p.z) + abs(m[2][3]);

        (
            self.transform_point(p),
            gamma(3) * vector3(x_abs_sum, y_abs_sum, z_abs_sum),
        )
    }

    /// Using the original point passed to `transform_point` and the error
    /// returned by `transform_point_error` returns the absolute error in the
    /// result.
    ///
    /// * `p`       - The point passed to `transform_point_with_error`.
    /// * `p_error` - The absolute error computed in `transform_point_with_error`.
    pub fn transform_point_abs_error(&self, p: &Point3f, p_error: &Vector3f) -> Vector3f {
        // Calculate the absolute error in the result based on its own error.
        let m = &self.m;
        let gamma_3 = gamma(3);

        let abs_error_x = (gamma_3 + 1.0)
            * (abs(m[0][0]) * p_error.x + abs(m[0][1]) * p_error.y + abs(m[0][2]) * p_error.z)
            + gamma_3
                * (abs(m[0][0] * p.x) + abs(m[0][1] * p.y) + abs(m[0][2] * p.z) + abs(m[0][3]));

        let abs_error_y = (gamma_3 + 1.0)
            * (abs(m[1][0]) * p_error.x + abs(m[1][1]) * p_error.y + abs(m[1][2]) * p_error.z)
            + gamma_3
                * (abs(m[1][0] * p.x) + abs(m[1][1] * p.y) + abs(m[1][2] * p.z) + abs(m[1][3]));

        let abs_error_z = (gamma_3 + 1.0)
            * (abs(m[2][0]) * p_error.x + abs(m[2][1]) * p_error.y + abs(m[2][2]) * p_error.z)
            + gamma_3
                * (abs(m[2][0] * p.x) + abs(m[2][1] * p.y) + abs(m[2][2] * p.z) + abs(m[2][3]));

        vector3(abs_error_x, abs_error_y, abs_error_z)
    }

    /// Applies transformation to a given vector.
    ///
    /// * `v` - The vector.
    pub fn transform_vector(&self, v: &Vector3f) -> Vector3f {
        let m = &self.m;
        vector3(
            m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z,
        )
    }

    /// Returns the transformed vector with absolute error due to applying the
    /// transformation to a vector.
    ///
    /// * `v` - The vector.
    pub fn transform_vector_with_error(&self, v: &Vector3f) -> (Vector3f, Vector3f) {
        // Compute absolute error for transformed vector
        let m = &self.m;

        let x_abs_err = abs(m[0][0] * v.x) + abs(m[0][1] * v.y) + abs(m[0][2] * v.z);
        let y_abs_err = abs(m[1][0] * v.x) + abs(m[1][1] * v.y) + abs(m[1][2] * v.z);
        let z_abs_err = abs(m[2][0] * v.x) + abs(m[2][1] * v.y) + abs(m[2][2] * v.z);

        (
            self.transform_vector(v),
            gamma(3) * vector3(x_abs_err, y_abs_err, z_abs_err),
        )
    }

    /// Using the original vector passed to `transform_vector` and the error
    /// returned by `transform_vector_error` returns the absolute error in the
    /// result.
    ///
    /// * `v`       - The vector passed to `transform_vector_with_error`.
    /// * `v_error` - The absolute error computed in `transform_vector_with_error`.
    pub fn transform_vector_abs_error(&self, v: &Vector3f, v_error: &Vector3f) -> Vector3f {
        // Calculate the absolute error in the result based on its own error.
        let m = &self.m;
        let gamma_3 = gamma(3);

        let abs_error_x = (gamma_3 + 1.0)
            * (abs(m[0][0]) * v_error.x + abs(m[0][1]) * v_error.y + abs(m[0][2]) * v_error.z)
            + gamma_3 * (abs(m[0][0] * v.x) + abs(m[0][1] * v.y) + abs(m[0][2] * v.z));

        let abs_error_y = (gamma_3 + 1.0)
            * (abs(m[1][0]) * v_error.x + abs(m[1][1]) * v_error.y + abs(m[1][2]) * v_error.z)
            + gamma_3 * (abs(m[1][0] * v.x) + abs(m[1][1] * v.y) + abs(m[1][2] * v.z));

        let abs_error_z = (gamma_3 + 1.0)
            * (abs(m[2][0]) * v_error.x + abs(m[2][1]) * v_error.y + abs(m[2][2]) * v_error.z)
            + gamma_3 * (abs(m[2][0] * v.x) + abs(m[2][1] * v.y) + abs(m[2][2] * v.z));

        vector3(abs_error_x, abs_error_y, abs_error_z)
    }

    /// Applies transformation to a given normal.
    ///
    /// * `n` - The normal.
    pub fn transform_normal(&self, n: &Normal3f) -> Normal3f {
        let m_inv = &self.m_inv.m;
        normal3(
            m_inv[0][0] * n.x + m_inv[1][0] * n.y + m_inv[2][0] * n.z,
            m_inv[0][1] * n.x + m_inv[1][1] * n.y + m_inv[2][1] * n.z,
            m_inv[0][2] * n.x + m_inv[1][2] * n.y + m_inv[2][2] * n.z,
        )
    }

    /// Applies transformation to a given ray.
    ///
    /// * `r` - The ray.
    pub fn transform_ray(&self, r: &Ray) -> Ray {
        let (mut o, o_error) = self.transform_point_with_error(&r.o);
        let d = self.transform_vector(&r.d);

        // Offset ray origin to edge of error bounds and compute t_max.
        let length_squared = d.length_squared();
        let mut t_max = r.t_max;
        if length_squared > 0.0 {
            let dt = d.abs().dot(&o_error) / length_squared;
            o += d * dt;
            t_max -= dt;
        }

        // Handle differentials.
        if let Some(diff) = r.differentials {
            let td = ray_differential(
                self.transform_point(&diff.rx_origin),
                self.transform_point(&diff.ry_origin),
                self.transform_vector(&diff.rx_direction),
                self.transform_vector(&diff.ry_direction),
            );
            ray_with_differentials(o, d, t_max, r.time, td, r.medium.clone())
        } else {
            ray(o, d, t_max, r.time, r.medium.clone())
        }
    }

    /// Returns the transformed ray with absolute errors due to applying the
    /// transformation to its origin and direction.
    ///
    /// * `r` - The ray.
    pub fn transform_ray_with_error(&self, r: &Ray) -> (Ray, Vector3f, Vector3f) {
        let (mut o, o_error) = self.transform_point_with_error(&r.o);
        let (d, d_error) = self.transform_vector_with_error(&r.d);

        // Offset ray origin to edge of error bounds.
        let length_sqquared = d.length_squared();
        if length_sqquared > 0.0 {
            let dt = d.abs().dot(&o_error) / length_sqquared;
            o += d * dt;
        }

        // Handle differentials.
        if let Some(diff) = r.differentials {
            let td = ray_differential(
                self.transform_point(&diff.rx_origin),
                self.transform_point(&diff.ry_origin),
                self.transform_vector(&diff.rx_direction),
                self.transform_vector(&diff.ry_direction),
            );

            (
                ray_with_differentials(o, d, r.t_max, r.time, td, r.medium.clone()),
                o_error,
                d_error,
            )
        } else {
            (
                ray(o, d, r.t_max, r.time, r.medium.clone()),
                o_error,
                d_error,
            )
        }
    }

    /// Transforms the ray taking into account absolute errors due to applying the
    /// transformation to its origin and direction.
    ///
    /// * `r` - The ray.
    pub fn transform_ray_with_abs_error(&self, r: &Ray) -> (Ray, Vector3f, Vector3f) {
        let (mut tr, o_error_in, d_error_in) = self.transform_ray_with_error(r);

        // Calculate error in result for origin and direction.
        let o_error_out = self.transform_point_abs_error(&tr.o, &o_error_in);
        let d_error_out = self.transform_vector_abs_error(&tr.d, &d_error_in);
        // Offset ray origin to edge of error bounds.
        let length_sqquared = tr.d.length_squared();
        if length_sqquared > 0.0 {
            let dt = tr.d.abs().dot(&o_error_out) / length_sqquared;
            tr.o += tr.d * dt;
        }

        (tr, o_error_out, d_error_out)
    }

    /// Applies transformation to a given bounding box.
    ///
    /// * `b` - The bounding box.
    pub fn transform_bounds(&self, b: &Bounds3f) -> Bounds3f {
        Bounds3::from(self.transform_point(&point3(b.p_min.x, b.p_min.y, b.p_min.z)))
            .union(&self.transform_point(&point3(b.p_max.x, b.p_min.y, b.p_min.z)))
            .union(&self.transform_point(&point3(b.p_min.x, b.p_max.y, b.p_min.z)))
            .union(&self.transform_point(&point3(b.p_min.x, b.p_min.y, b.p_max.z)))
            .union(&self.transform_point(&point3(b.p_min.x, b.p_max.y, b.p_max.z)))
            .union(&self.transform_point(&point3(b.p_max.x, b.p_max.y, b.p_min.z)))
            .union(&self.transform_point(&point3(b.p_max.x, b.p_min.y, b.p_max.z)))
            .union(&self.transform_point(&point3(b.p_max.x, b.p_max.y, b.p_max.z)))
    }

    /// Applies transformation to a given surface interaction.
    ///
    /// * `si` - The surface interaction.
    pub fn transform_surface_interaction(&self, si: &SurfaceInteraction) -> SurfaceInteraction {
        // Transform p and p_error in SurfaceInteraction
        let (p, p_error) = self.transform_point_with_error(&si.hit.p);

        // Transform remaining members of SurfaceInteraction
        let mut si = surface_interaction(
            p,
            p_error,
            si.uv,
            self.transform_vector(&si.hit.wo).normalize(),
            self.transform_vector(&si.dpdu),
            self.transform_vector(&si.dpdv),
            self.transform_normal(&si.dndu),
            self.transform_normal(&si.dndv),
            si.hit.time,
            si.shape.clone(),
        );

        // Transform n in SurfaceInteraction.hit
        let n = self.transform_normal(&si.hit.n).normalize();
        si.hit.n = n;

        // Handle transformations for shading parameters..
        si.shading = shading(
            self.transform_normal(&si.shading.n).normalize(),
            self.transform_vector(&si.shading.dpdu),
            self.transform_vector(&si.shading.dpdv),
            self.transform_normal(&si.shading.dndu),
            self.transform_normal(&si.shading.dndv),
        );
        si.shading.n = si.shading.n.face_forward(&Vector3::from(n));

        si
    }

    /// Returns `true` if the transformation changes the handedness of the
    /// coordinate system.
    pub fn swaps_handedness(&self) -> bool {
        let m = &self.m;
        let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        det < 0.0
    }
}

impl From<Matrix4x4> for Transform {
    /// Create a transformation from a 4x4 matrix.
    ///
    /// * `m` - A matrix representing a transformation.
    fn from(m: Matrix4x4) -> Self {
        Transform {
            m,
            m_inv: m.inverse(),
        }
    }
}

impl PartialOrd for Transform {
    fn partial_cmp(&self, other: &Transform) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Transform {
    fn cmp(&self, other: &Transform) -> Ordering {
        for i in 0..4 {
            for j in 0..4 {
                if self.m[i][j] < other.m[i][j] {
                    return Ordering::Less;
                }
                if self.m[i][j] > other.m[i][j] {
                    return Ordering::Greater;
                }
            }
        }
        Ordering::Equal
    }
}

impl PartialEq for Transform {
    fn eq(&self, other: &Self) -> bool {
        self.m == other.m
    }
}

impl Eq for Transform {}

impl Mul<Transform> for Transform {
    type Output = Self;

    /// Composes this transformation with another one. The resulting transform
    /// is the same as applying `self` then `rhs`.
    ///
    /// * `rhs` - The transformation to compose.
    fn mul(self, rhs: Self) -> Self {
        Transform {
            m: self.m * rhs.m,
            m_inv: rhs.m_inv * self.m_inv,
        }
    }
}
