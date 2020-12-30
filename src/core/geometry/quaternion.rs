//! Quaternions

#![allow(dead_code)]
use super::common::*;
use crate::core::geometry::*;
use crate::core::pbrt::*;
use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// A quaternion
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Quaternion {
    /// The `x`, `y`, `z` components represented as a vector.
    pub v: Vector3f,

    /// The real component `w`.
    pub w: Float,
}

impl Quaternion {
    /// Create a new quaternion.
    ///
    /// * `v` - The `x`, `y`, `z` components represented as a vector.
    /// * `w` - The real component `w`.
    pub fn new(v: Vector3f, w: Float) -> Self {
        Self { v, w }
    }

    /// Normlizes the quaternion by dividing each component by its length.
    pub fn normalize(&self) -> Self {
        *self / self.length()
    }

    /// Returns the square of the length of the quaternion which is the
    /// inner product with itself.
    pub fn length_squared(&self) -> Float {
        self.dot(self)
    }

    /// Returns the length of the quaternion which is square root of the inner
    /// product with itself.
    pub fn length(&self) -> Float {
        self.length_squared().sqrt()
    }

    /// Interpolate between this and another quaternion using sqpherical linear
    /// interpolation.
    ///
    /// * `t` - The interpolation parameter.
    /// * `q` - The other quaternion.
    pub fn slerp(&self, t: Float, q: Self) -> Self {
        let cos_theta = self.dot(&q);
        if cos_theta > 0.9995 {
            // Quaternions are nearly parallel. Use linear interpolation to
            // avoid numerical instability.
            ((1.0 - t) * *self + t * q).normalize()
        } else {
            //  Compute the orthogonal quaternion `qperp`.
            let theta = clamp(cos_theta, -1.0, 1.0).acos();
            let thetap = theta * t;
            let qperp = (q - *self * cos_theta).normalize();

            // Compute the interpolated quaternion.
            *self * thetap.cos() + qperp * thetap.sin()
        }
    }
}

impl Default for Quaternion {
    /// Returns the default quaternion [0, 0, 0, 1].
    fn default() -> Self {
        Self {
            v: Vector3f::zero(),
            w: 1.0,
        }
    }
}

impl fmt::Display for Quaternion {
    /// Formats the value using the given formatter.
    ///
    /// * `f` - Formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}, {}, {}, {}]", self.v.x, self.v.y, self.v.z, self.w)
    }
}

impl From<Transform> for Quaternion {
    /// Returns a quaternion representing a rotational transform.
    ///
    /// * `t` - The transform.
    fn from(t: Transform) -> Self {
        let m = &t.m;
        let trace = m[0][0] + m[1][1] + m[2][2];
        if trace > 0.0 {
            // Compute w from matrix trace, then xyz
            // 4w^2 = m[0][0] + m[1][1] + m[2][2] + m[3][3] (but m[3][3] == 1)
            let mut s = (trace + 1.0).sqrt();

            let w = s / 2.0;

            s = 0.5 / s;

            let v = Vector3f::new(
                (m[2][1] - m[1][2]) * s,
                (m[0][2] - m[2][0]) * s,
                (m[1][0] - m[0][1]) * s,
            );

            Self::new(v, w)
        } else {
            // Compute largest of x, y, or z, then remaining components
            let nxt = [1, 2, 0];
            let mut q = [0.0; 3];

            let mut i = 0;
            if m[1][1] > m[0][0] {
                i = 1;
            }
            if m[2][2] > m[i][i] {
                i = 2;
            }

            let j = nxt[i];
            let k = nxt[j];

            let mut s = ((m[i][i] - (m[j][j] + m[k][k])) + 1.0).sqrt();
            q[i] = s * 0.5;

            if s != 0.0 {
                s = 0.5 / s;
            }

            let w = (m[k][j] - m[j][k]) * s;

            q[j] = (m[j][i] + m[i][j]) * s;
            q[k] = (m[k][i] + m[i][k]) * s;

            let v = Vector3f::new(q[0], q[1], q[2]);

            Self::new(v, w)
        }
    }
}

impl From<Quaternion> for Transform {
    /// Returns a rotation transform from a quaternion.
    ///
    /// * `q` - The quaternion.
    fn from(q: Quaternion) -> Transform {
        let xx = q.v.x * q.v.x;
        let yy = q.v.y * q.v.y;
        let zz = q.v.z * q.v.z;
        let xy = q.v.x * q.v.y;
        let xz = q.v.x * q.v.z;
        let yz = q.v.y * q.v.z;
        let wx = q.v.x * q.w;
        let wy = q.v.y * q.w;
        let wz = q.v.z * q.w;

        let mut m = IDENTITY_MATRIX;
        m.m[0][0] = 1.0 - 2.0 * (yy + zz);
        m.m[0][1] = 2.0 * (xy + wz);
        m.m[0][2] = 2.0 * (xz - wy);
        m.m[1][0] = 2.0 * (xy - wz);
        m.m[1][1] = 1.0 - 2.0 * (xx + zz);
        m.m[1][2] = 2.0 * (yz + wx);
        m.m[2][0] = 2.0 * (xz + wy);
        m.m[2][1] = 2.0 * (yz - wx);
        m.m[2][2] = 1.0 - 2.0 * (xx + yy);

        // Transpose since we are left-handed.
        Transform {
            m: m.transpose(),
            m_inv: m,
        }
    }
}

impl Add<Quaternion> for Quaternion {
    type Output = Self;

    /// Adds the components of the given quaternion.
    ///
    /// * `other` - The quaternion to add.
    fn add(self, other: Self) -> Self::Output {
        Self::Output::new(self.v + other.v, self.w + other.w)
    }
}

impl AddAssign<Quaternion> for Quaternion {
    /// Performs the `+=` operation.
    ///
    /// * `other` - The quaternion to add.
    fn add_assign(&mut self, other: Self) {
        *self = Self::new(self.v + other.v, self.w + other.w)
    }
}

impl Sub<Quaternion> for Quaternion {
    type Output = Self;

    /// Subtracts the components of the given quaternion.
    ///
    /// * `other` - The quaternion to subtract.
    fn sub(self, other: Self) -> Self::Output {
        Self::new(self.v - other.v, self.w - other.w)
    }
}

impl SubAssign<Quaternion> for Quaternion {
    /// Performs the `-=` operation.
    ///
    /// * `other` - The quaternion to subtract.
    fn sub_assign(&mut self, other: Self) {
        *self = Self::new(self.v - other.v, self.w - other.w)
    }
}

impl Mul<Float> for Quaternion {
    type Output = Self;

    /// Scales the components of the given quaternion.
    ///
    /// * `f` - The scaling factor.
    fn mul(self, f: Float) -> Self::Output {
        Self::Output::new(f * self.v, f * self.w)
    }
}

impl Mul<Quaternion> for Float {
    type Output = Quaternion;

    /// Scales the components of the given quaternion.
    ///
    /// * `q` - The quaternion to scale.
    fn mul(self, q: Quaternion) -> Self::Output {
        Self::Output::new(self * q.v, self * q.w)
    }
}

impl MulAssign<Float> for Quaternion {
    /// Performs the `*=` operation.
    ///
    /// * `f` - The scaling factor.
    fn mul_assign(&mut self, f: Float) {
        *self = Self::new(f * self.v, f * self.w)
    }
}

impl Div<Float> for Quaternion {
    type Output = Self;

    /// Scales the components of the given quaternion by 1/f.
    ///
    /// * `f` - The scaling factor.
    fn div(self, f: Float) -> Self::Output {
        Self::Output::new(self.v / f, self.w / f)
    }
}

impl DivAssign<Float> for Quaternion {
    /// Performs the `*=` operation.
    ///
    /// * `f` - The scaling factor.
    fn div_assign(&mut self, f: Float) {
        *self = Self::new(self.v / f, self.w / f)
    }
}

impl Neg for Quaternion {
    type Output = Self;

    /// Scales the components by -1.
    fn neg(self) -> Self::Output {
        Self::Output::new(-self.v, -self.w)
    }
}

impl Dot<Quaternion> for Quaternion {
    type Output = Float;

    /// Returns the inner product with another quaternion.
    ///
    /// * `other` - The other quaternion.
    fn dot(&self, other: &Quaternion) -> Float {
        self.v.dot(&other.v) + self.w * other.w
    }
}
