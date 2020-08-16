//! 4x4 Matrix

#![allow(dead_code)]
use super::{abs, Float};
use std::ops::{Index, Mul};

/// A 4x4 vector containing Float values.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Matrix4x4 {
    /// Stores a 2-D array of Float
    pub m: [[Float; 4]; 4],
}

/// Zero matrix.
pub const ZERO_MATRIX: Matrix4x4 = Matrix4x4 { m: [[0.0; 4]; 4] };

/// Ientity matrix.
pub const IDENTITY_MATRIX: Matrix4x4 = Matrix4x4 {
    m: [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ],
};

/// Create a 4x4 matrix using the following order of the parameters:
///
/// * `t00`, `t01`, `t02`, `t03` - Row 1
/// * `t10`, `t11`, `t12`, `t13` - Row 2
/// * `t20`, `t21`, `t22`, `t23` - Row 3
/// * `t30`, `t31`, `t32`, `t33` - Row 4
#[rustfmt::skip]
pub fn matrix4x4(
    t00: Float, t01: Float, t02: Float, t03: Float,
    t10: Float, t11: Float, t12: Float, t13: Float,
    t20: Float, t21: Float, t22: Float, t23: Float,
    t30: Float, t31: Float, t32: Float, t33: Float,
) -> Matrix4x4 {
    Matrix4x4 {
        m: [
            [t00, t01, t02, t03],
            [t10, t11, t12, t13],
            [t20, t21, t22, t23],
            [t30, t31, t32, t33],
        ],
    }
}

#[rustfmt::skip]
impl Matrix4x4 {
    /// Returns the transpose of the matrix.
    pub fn transpose(&self) -> Matrix4x4 {
        matrix4x4(
            self.m[0][0], self.m[1][0], self.m[2][0], self.m[3][0],
            self.m[0][1], self.m[1][1], self.m[2][1], self.m[3][1],
            self.m[0][2], self.m[1][2], self.m[2][2], self.m[3][2],
            self.m[0][3], self.m[1][3], self.m[2][3], self.m[3][3],
        )
    }

    /// Returns the inverse of the matrix using numerically stable Gauss-Jordan
    /// elimination.
    ///
    /// The function will panic if the matrix is singular.
    pub fn inverse(&self) -> Matrix4x4 {
        let mut indxc = [0; 4];
        let mut indxr = [0; 4];
        let mut ipiv = [0; 4];

        let mut minv: [[Float; 4]; 4] = [[0.0; 4]; 4];
        minv.copy_from_slice(&self.m);

        for i in 0..4 {
            let mut irow = 0;
            let mut icol = 0;
            let mut big: Float = 0.0;

            // Choose pivot
            for j in 0..4 {
                if ipiv[j] != 1 {
                    for k in 0..4 {
                        if ipiv[k] == 0 {
                            let abs_minv = abs(minv[j][k]);
                            if abs_minv >= big {
                                big = abs_minv;
                                irow = j;
                                icol = k;
                            }
                        } else if ipiv[k] > 1 {
                            panic!("Singular matrix in MatrixInvert");
                        }
                    }
                }
            }

            ipiv[icol] += 1;

            // Swap rows _irow_ and _icol_ for pivot
            if irow != icol {
                for k in 0..4 {
                    let tmp = minv[irow][k];
                    minv[irow][k] = minv[icol][k];
                    minv[icol][k] = tmp;
                }
            }

            indxr[i] = irow;
            indxc[i] = icol;
            if minv[icol][icol] == 0.0 {
                panic!("Singular matrix in MatrixInvert");
            }

            // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
            let pivinv = 1.0 / minv[icol][icol];
            minv[icol][icol] = 1.0;
            for j in 0..4 {
                minv[icol][j] *= pivinv;
            }

            // Subtract this row from others to zero out their columns
            for j in 0..4 {
                if j != icol {
                    let save = minv[j][icol];
                    minv[j][icol] = 0.0;
                    for k in 0..4 {
                        minv[j][k] -= minv[icol][k] * save;
                    }
                }
            }
        }
        // Swap columns to reflect permutation
        for j in (0..4).rev() {
            if indxr[j] != indxc[j] {
                for k in 0..4 {
                    let tmp = minv[k][indxr[j]];
                    minv[k][indxr[j]] = minv[k][indxc[j]];
                    minv[k][indxc[j]] = tmp;
                }
            }
        }

        Matrix4x4 { m: minv }
    }
}

impl Default for Matrix4x4 {
    /// Returns the default as identity matrix.
    fn default() -> Self {
        IDENTITY_MATRIX
    }
}

impl Mul<Matrix4x4> for Matrix4x4 {
    type Output = Matrix4x4;

    /// Post-multiply the given matrix.
    ///
    /// * `other` - The other matrix
    fn mul(self, other: Matrix4x4) -> Self::Output {
        let mut m = Matrix4x4::default();

        for i in 0..4 {
            for j in 0..4 {
                m.m[i][j] = self.m[i][0] * other.m[0][j]
                    + self.m[i][1] * other.m[1][j]
                    + self.m[i][2] * other.m[2][j]
                    + self.m[i][3] * other.m[3][j];
            }
        }

        m
    }
}

impl Index<usize> for Matrix4x4 {
    type Output = [Float; 4];

    /// Index the matrix row. The column can be further indexed from the
    /// returned result.
    ///
    /// * `row` - Row
    fn index(&self, row: usize) -> &Self::Output {
        assert!(row < 4, "matrix row not in [0, 3]");
        &self.m[row]
    }
}

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

#[cfg(test)]
#[macro_use]
mod tests {
    use super::*;
    use float_cmp::*;
    use proptest::prelude::*;

    #[test]
    #[should_panic]
    fn inverse_panics_when_matrix_is_zero() {
        let _ = ZERO_MATRIX.inverse();
    }

    #[test]
    fn inverse_returns_identity_when_matrix_is_idenitity() {
        assert_eq!(IDENTITY_MATRIX.inverse(), IDENTITY_MATRIX);
    }

    proptest! {
        #[test]
        #[should_panic]
        fn inverse_panics_when_matrix_is_singular(
            a in 0.0..10.0f32, b in 0.0..10.0f32, c in 0.0..10.0f32,
        ) {
            let  _ = Matrix4x4 {
                m: [
                    [1.0, 0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0],
                    [  a,   b,   c, 0.0],
                ],
            }.inverse();
        }

        #[test]
        fn inverse_returns_matrix_when_matrix_is_non_singular(
            a in 0.001..10.0f32, b in 0.001..10.0f32, c in 0.001..10.0f32, d in 0.001..10.0f32,
        ) {
            let mat = Matrix4x4 {
                m: [
                    [  a, 0.0, 0.0, 0.0],
                    [0.0,   b, 0.0, 0.0],
                    [0.0, 0.0,   c, 0.0],
                    [0.0, 0.0, 0.0,   d],
                ],
            };

            let prod = mat * mat.inverse();
            for j in 0..4 {
                for i in 0..4 {
                    prop_assert!(approx_eq!(
                            Float,
                            prod.m[i][j],
                            IDENTITY_MATRIX.m[i][j],
                            epsilon = 0.0001
                    ));
                }
            }

            let prod = mat.inverse() * mat;
            for j in 0..4 {
                for i in 0..4 {
                    prop_assert!(approx_eq!(
                            Float,
                            prod.m[i][j],
                            IDENTITY_MATRIX.m[i][j],
                            epsilon = 0.0001
                    ));
                }
            }
        }
    }
}
