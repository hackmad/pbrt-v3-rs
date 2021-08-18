//! Animated Transformations

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::pbrt::*;
use std::sync::Arc;

/// AnimatedTransform implements keyframe transformation interpolation.
#[derive(Clone, Debug)]
pub struct AnimatedTransform {
    /// Transform representing starting point in animation.
    pub start_transform: ArcTransform,

    /// Transform representing ending point in animation.
    pub end_transform: ArcTransform,

    /// Start time.
    pub start_time: Float,

    /// End time.
    pub end_time: Float,

    /// Indicates if `start_transform` == `end_transform`.
    actually_animated: bool,

    /// The translation component of decomposed transformation matrix.
    t: [Vector3f; 2],

    /// The rotation component of decomposed transformation matrix.
    r: [Quaternion; 2],

    /// The scaling component of decomposed transformation matrix.
    s: [Matrix4x4; 2],

    /// Indicates the animated transformation has rotation.
    has_rotation: bool,

    /// Stores derivative information to bind the motion of a point.
    c1: [DerivativeTerm; 3],
    c2: [DerivativeTerm; 3],
    c3: [DerivativeTerm; 3],
    c4: [DerivativeTerm; 3],
    c5: [DerivativeTerm; 3],
}

impl AnimatedTransform {
    /// Create an animated transformation.
    ///
    /// * `start_transform` - Transform representing starting point in animation.
    /// * `end_transform`   - transform representing ending point in animation.
    /// * `start_time`      - start time.
    /// * `end_time`        - End time.
    pub fn new(
        start_transform: ArcTransform,
        end_transform: ArcTransform,
        start_time: Float,
        end_time: Float,
    ) -> Self {
        let start_transform = Arc::clone(&start_transform);
        let end_transform = Arc::clone(&end_transform);
        let actually_animated = start_transform != end_transform;
        let mut t = [Vector3f::default(); 2];
        let mut r = [Quaternion::default(); 2];
        let mut s = [Matrix4x4::default(); 2];
        let mut has_rotation = false;
        let mut c1 = [DerivativeTerm::default(); 3];
        let mut c2 = [DerivativeTerm::default(); 3];
        let mut c3 = [DerivativeTerm::default(); 3];
        let mut c4 = [DerivativeTerm::default(); 3];
        let mut c5 = [DerivativeTerm::default(); 3];

        if actually_animated {
            decompose(&start_transform.m, &mut t[0], &mut r[0], &mut s[0]);
            decompose(&end_transform.m, &mut t[1], &mut r[1], &mut s[1]);

            // Flip R[1] if needed to select shortest path
            if r[0].dot(&r[1]) < 0.0 {
                r[1] = -r[1];
            }

            let cos_theta = r[0].dot(&r[1]);
            has_rotation = cos_theta < 0.9995;

            // Compute terms of motion derivative function
            if has_rotation {
                let theta = clamp(cos_theta, -1.0, 1.0).acos();
                let qperp = (r[1] - r[0] * cos_theta).normalize();

                let t0x = t[0].x;
                let t0y = t[0].y;
                let t0z = t[0].z;
                let t1x = t[1].x;
                let t1y = t[1].y;
                let t1z = t[1].z;
                let q1x = r[0].v.x;
                let q1y = r[0].v.y;
                let q1z = r[0].v.z;
                let q1w = r[0].w;
                let qperpx = qperp.v.x;
                let qperpy = qperp.v.y;
                let qperpz = qperp.v.z;
                let qperpw = qperp.w;
                let s000 = s[0].m[0][0];
                let s001 = s[0].m[0][1];
                let s002 = s[0].m[0][2];
                let s010 = s[0].m[1][0];
                let s011 = s[0].m[1][1];
                let s012 = s[0].m[1][2];
                let s020 = s[0].m[2][0];
                let s021 = s[0].m[2][1];
                let s022 = s[0].m[2][2];
                let s100 = s[1].m[0][0];
                let s101 = s[1].m[0][1];
                let s102 = s[1].m[0][2];
                let s110 = s[1].m[1][0];
                let s111 = s[1].m[1][1];
                let s112 = s[1].m[1][2];
                let s120 = s[1].m[2][0];
                let s121 = s[1].m[2][1];
                let s122 = s[1].m[2][2];

                c1[0] = derivative_term(
                    -t0x + t1x,
                    (-1.0 + q1y * q1y + q1z * q1z + qperpy * qperpy + qperpz * qperpz) * s000
                        + q1w * q1z * s010
                        - qperpx * qperpy * s010
                        + qperpw * qperpz * s010
                        - q1w * q1y * s020
                        - qperpw * qperpy * s020
                        - qperpx * qperpz * s020
                        + s100
                        - q1y * q1y * s100
                        - q1z * q1z * s100
                        - qperpy * qperpy * s100
                        - qperpz * qperpz * s100
                        - q1w * q1z * s110
                        + qperpx * qperpy * s110
                        - qperpw * qperpz * s110
                        + q1w * q1y * s120
                        + qperpw * qperpy * s120
                        + qperpx * qperpz * s120
                        + q1x * (-(q1y * s010) - q1z * s020 + q1y * s110 + q1z * s120),
                    (-1.0 + q1y * q1y + q1z * q1z + qperpy * qperpy + qperpz * qperpz) * s001
                        + q1w * q1z * s011
                        - qperpx * qperpy * s011
                        + qperpw * qperpz * s011
                        - q1w * q1y * s021
                        - qperpw * qperpy * s021
                        - qperpx * qperpz * s021
                        + s101
                        - q1y * q1y * s101
                        - q1z * q1z * s101
                        - qperpy * qperpy * s101
                        - qperpz * qperpz * s101
                        - q1w * q1z * s111
                        + qperpx * qperpy * s111
                        - qperpw * qperpz * s111
                        + q1w * q1y * s121
                        + qperpw * qperpy * s121
                        + qperpx * qperpz * s121
                        + q1x * (-(q1y * s011) - q1z * s021 + q1y * s111 + q1z * s121),
                    (-1.0 + q1y * q1y + q1z * q1z + qperpy * qperpy + qperpz * qperpz) * s002
                        + q1w * q1z * s012
                        - qperpx * qperpy * s012
                        + qperpw * qperpz * s012
                        - q1w * q1y * s022
                        - qperpw * qperpy * s022
                        - qperpx * qperpz * s022
                        + s102
                        - q1y * q1y * s102
                        - q1z * q1z * s102
                        - qperpy * qperpy * s102
                        - qperpz * qperpz * s102
                        - q1w * q1z * s112
                        + qperpx * qperpy * s112
                        - qperpw * qperpz * s112
                        + q1w * q1y * s122
                        + qperpw * qperpy * s122
                        + qperpx * qperpz * s122
                        + q1x * (-(q1y * s012) - q1z * s022 + q1y * s112 + q1z * s122),
                );

                c2[0] = derivative_term(
                    0.0,
                    -(qperpy * qperpy * s000) - qperpz * qperpz * s000 + qperpx * qperpy * s010
                        - qperpw * qperpz * s010
                        + qperpw * qperpy * s020
                        + qperpx * qperpz * s020
                        + q1y * q1y * (s000 - s100)
                        + q1z * q1z * (s000 - s100)
                        + qperpy * qperpy * s100
                        + qperpz * qperpz * s100
                        - qperpx * qperpy * s110
                        + qperpw * qperpz * s110
                        - qperpw * qperpy * s120
                        - qperpx * qperpz * s120
                        + 2.0 * q1x * qperpy * s010 * theta
                        - 2.0 * q1w * qperpz * s010 * theta
                        + 2.0 * q1w * qperpy * s020 * theta
                        + 2.0 * q1x * qperpz * s020 * theta
                        + q1y
                            * (q1x * (-s010 + s110)
                                + q1w * (-s020 + s120)
                                + 2.0
                                    * (-2.0 * qperpy * s000 + qperpx * s010 + qperpw * s020)
                                    * theta)
                        + q1z
                            * (q1w * (s010 - s110) + q1x * (-s020 + s120)
                                - 2.0
                                    * (2.0 * qperpz * s000 + qperpw * s010 - qperpx * s020)
                                    * theta),
                    -(qperpy * qperpy * s001) - qperpz * qperpz * s001 + qperpx * qperpy * s011
                        - qperpw * qperpz * s011
                        + qperpw * qperpy * s021
                        + qperpx * qperpz * s021
                        + q1y * q1y * (s001 - s101)
                        + q1z * q1z * (s001 - s101)
                        + qperpy * qperpy * s101
                        + qperpz * qperpz * s101
                        - qperpx * qperpy * s111
                        + qperpw * qperpz * s111
                        - qperpw * qperpy * s121
                        - qperpx * qperpz * s121
                        + 2.0 * q1x * qperpy * s011 * theta
                        - 2.0 * q1w * qperpz * s011 * theta
                        + 2.0 * q1w * qperpy * s021 * theta
                        + 2.0 * q1x * qperpz * s021 * theta
                        + q1y
                            * (q1x * (-s011 + s111)
                                + q1w * (-s021 + s121)
                                + 2.0
                                    * (-2.0 * qperpy * s001 + qperpx * s011 + qperpw * s021)
                                    * theta)
                        + q1z
                            * (q1w * (s011 - s111) + q1x * (-s021 + s121)
                                - 2.0
                                    * (2.0 * qperpz * s001 + qperpw * s011 - qperpx * s021)
                                    * theta),
                    -(qperpy * qperpy * s002) - qperpz * qperpz * s002 + qperpx * qperpy * s012
                        - qperpw * qperpz * s012
                        + qperpw * qperpy * s022
                        + qperpx * qperpz * s022
                        + q1y * q1y * (s002 - s102)
                        + q1z * q1z * (s002 - s102)
                        + qperpy * qperpy * s102
                        + qperpz * qperpz * s102
                        - qperpx * qperpy * s112
                        + qperpw * qperpz * s112
                        - qperpw * qperpy * s122
                        - qperpx * qperpz * s122
                        + 2.0 * q1x * qperpy * s012 * theta
                        - 2.0 * q1w * qperpz * s012 * theta
                        + 2.0 * q1w * qperpy * s022 * theta
                        + 2.0 * q1x * qperpz * s022 * theta
                        + q1y
                            * (q1x * (-s012 + s112)
                                + q1w * (-s022 + s122)
                                + 2.0
                                    * (-2.0 * qperpy * s002 + qperpx * s012 + qperpw * s022)
                                    * theta)
                        + q1z
                            * (q1w * (s012 - s112) + q1x * (-s022 + s122)
                                - 2.0
                                    * (2.0 * qperpz * s002 + qperpw * s012 - qperpx * s022)
                                    * theta),
                );

                c3[0] = derivative_term(
                    0.0,
                    -2.0 * (q1x * qperpy * s010 - q1w * qperpz * s010
                        + q1w * qperpy * s020
                        + q1x * qperpz * s020
                        - q1x * qperpy * s110
                        + q1w * qperpz * s110
                        - q1w * qperpy * s120
                        - q1x * qperpz * s120
                        + q1y
                            * (-2.0 * qperpy * s000
                                + qperpx * s010
                                + qperpw * s020
                                + 2.0 * qperpy * s100
                                - qperpx * s110
                                - qperpw * s120)
                        + q1z
                            * (-2.0 * qperpz * s000 - qperpw * s010
                                + qperpx * s020
                                + 2.0 * qperpz * s100
                                + qperpw * s110
                                - qperpx * s120))
                        * theta,
                    -2.0 * (q1x * qperpy * s011 - q1w * qperpz * s011
                        + q1w * qperpy * s021
                        + q1x * qperpz * s021
                        - q1x * qperpy * s111
                        + q1w * qperpz * s111
                        - q1w * qperpy * s121
                        - q1x * qperpz * s121
                        + q1y
                            * (-2.0 * qperpy * s001
                                + qperpx * s011
                                + qperpw * s021
                                + 2.0 * qperpy * s101
                                - qperpx * s111
                                - qperpw * s121)
                        + q1z
                            * (-2.0 * qperpz * s001 - qperpw * s011
                                + qperpx * s021
                                + 2.0 * qperpz * s101
                                + qperpw * s111
                                - qperpx * s121))
                        * theta,
                    -2.0 * (q1x * qperpy * s012 - q1w * qperpz * s012
                        + q1w * qperpy * s022
                        + q1x * qperpz * s022
                        - q1x * qperpy * s112
                        + q1w * qperpz * s112
                        - q1w * qperpy * s122
                        - q1x * qperpz * s122
                        + q1y
                            * (-2.0 * qperpy * s002
                                + qperpx * s012
                                + qperpw * s022
                                + 2.0 * qperpy * s102
                                - qperpx * s112
                                - qperpw * s122)
                        + q1z
                            * (-2.0 * qperpz * s002 - qperpw * s012
                                + qperpx * s022
                                + 2.0 * qperpz * s102
                                + qperpw * s112
                                - qperpx * s122))
                        * theta,
                );

                c4[0] = derivative_term(
                    0.0,
                    -(q1x * qperpy * s010) + q1w * qperpz * s010
                        - q1w * qperpy * s020
                        - q1x * qperpz * s020
                        + q1x * qperpy * s110
                        - q1w * qperpz * s110
                        + q1w * qperpy * s120
                        + q1x * qperpz * s120
                        + 2.0 * q1y * q1y * s000 * theta
                        + 2.0 * q1z * q1z * s000 * theta
                        - 2.0 * qperpy * qperpy * s000 * theta
                        - 2.0 * qperpz * qperpz * s000 * theta
                        + 2.0 * qperpx * qperpy * s010 * theta
                        - 2.0 * qperpw * qperpz * s010 * theta
                        + 2.0 * qperpw * qperpy * s020 * theta
                        + 2.0 * qperpx * qperpz * s020 * theta
                        + q1y
                            * (-(qperpx * s010) - qperpw * s020
                                + 2.0 * qperpy * (s000 - s100)
                                + qperpx * s110
                                + qperpw * s120
                                - 2.0 * q1x * s010 * theta
                                - 2.0 * q1w * s020 * theta)
                        + q1z
                            * (2.0 * qperpz * s000 + qperpw * s010
                                - qperpx * s020
                                - 2.0 * qperpz * s100
                                - qperpw * s110
                                + qperpx * s120
                                + 2.0 * q1w * s010 * theta
                                - 2.0 * q1x * s020 * theta),
                    -(q1x * qperpy * s011) + q1w * qperpz * s011
                        - q1w * qperpy * s021
                        - q1x * qperpz * s021
                        + q1x * qperpy * s111
                        - q1w * qperpz * s111
                        + q1w * qperpy * s121
                        + q1x * qperpz * s121
                        + 2.0 * q1y * q1y * s001 * theta
                        + 2.0 * q1z * q1z * s001 * theta
                        - 2.0 * qperpy * qperpy * s001 * theta
                        - 2.0 * qperpz * qperpz * s001 * theta
                        + 2.0 * qperpx * qperpy * s011 * theta
                        - 2.0 * qperpw * qperpz * s011 * theta
                        + 2.0 * qperpw * qperpy * s021 * theta
                        + 2.0 * qperpx * qperpz * s021 * theta
                        + q1y
                            * (-(qperpx * s011) - qperpw * s021
                                + 2.0 * qperpy * (s001 - s101)
                                + qperpx * s111
                                + qperpw * s121
                                - 2.0 * q1x * s011 * theta
                                - 2.0 * q1w * s021 * theta)
                        + q1z
                            * (2.0 * qperpz * s001 + qperpw * s011
                                - qperpx * s021
                                - 2.0 * qperpz * s101
                                - qperpw * s111
                                + qperpx * s121
                                + 2.0 * q1w * s011 * theta
                                - 2.0 * q1x * s021 * theta),
                    -(q1x * qperpy * s012) + q1w * qperpz * s012
                        - q1w * qperpy * s022
                        - q1x * qperpz * s022
                        + q1x * qperpy * s112
                        - q1w * qperpz * s112
                        + q1w * qperpy * s122
                        + q1x * qperpz * s122
                        + 2.0 * q1y * q1y * s002 * theta
                        + 2.0 * q1z * q1z * s002 * theta
                        - 2.0 * qperpy * qperpy * s002 * theta
                        - 2.0 * qperpz * qperpz * s002 * theta
                        + 2.0 * qperpx * qperpy * s012 * theta
                        - 2.0 * qperpw * qperpz * s012 * theta
                        + 2.0 * qperpw * qperpy * s022 * theta
                        + 2.0 * qperpx * qperpz * s022 * theta
                        + q1y
                            * (-(qperpx * s012) - qperpw * s022
                                + 2.0 * qperpy * (s002 - s102)
                                + qperpx * s112
                                + qperpw * s122
                                - 2.0 * q1x * s012 * theta
                                - 2.0 * q1w * s022 * theta)
                        + q1z
                            * (2.0 * qperpz * s002 + qperpw * s012
                                - qperpx * s022
                                - 2.0 * qperpz * s102
                                - qperpw * s112
                                + qperpx * s122
                                + 2.0 * q1w * s012 * theta
                                - 2.0 * q1x * s022 * theta),
                );

                c5[0] = derivative_term(
                    0.0,
                    2.0 * (qperpy * qperpy * s000 + qperpz * qperpz * s000
                        - qperpx * qperpy * s010
                        + qperpw * qperpz * s010
                        - qperpw * qperpy * s020
                        - qperpx * qperpz * s020
                        - qperpy * qperpy * s100
                        - qperpz * qperpz * s100
                        + q1y * q1y * (-s000 + s100)
                        + q1z * q1z * (-s000 + s100)
                        + qperpx * qperpy * s110
                        - qperpw * qperpz * s110
                        + q1y * (q1x * (s010 - s110) + q1w * (s020 - s120))
                        + qperpw * qperpy * s120
                        + qperpx * qperpz * s120
                        + q1z * (-(q1w * s010) + q1x * s020 + q1w * s110 - q1x * s120))
                        * theta,
                    2.0 * (qperpy * qperpy * s001 + qperpz * qperpz * s001
                        - qperpx * qperpy * s011
                        + qperpw * qperpz * s011
                        - qperpw * qperpy * s021
                        - qperpx * qperpz * s021
                        - qperpy * qperpy * s101
                        - qperpz * qperpz * s101
                        + q1y * q1y * (-s001 + s101)
                        + q1z * q1z * (-s001 + s101)
                        + qperpx * qperpy * s111
                        - qperpw * qperpz * s111
                        + q1y * (q1x * (s011 - s111) + q1w * (s021 - s121))
                        + qperpw * qperpy * s121
                        + qperpx * qperpz * s121
                        + q1z * (-(q1w * s011) + q1x * s021 + q1w * s111 - q1x * s121))
                        * theta,
                    2.0 * (qperpy * qperpy * s002 + qperpz * qperpz * s002
                        - qperpx * qperpy * s012
                        + qperpw * qperpz * s012
                        - qperpw * qperpy * s022
                        - qperpx * qperpz * s022
                        - qperpy * qperpy * s102
                        - qperpz * qperpz * s102
                        + q1y * q1y * (-s002 + s102)
                        + q1z * q1z * (-s002 + s102)
                        + qperpx * qperpy * s112
                        - qperpw * qperpz * s112
                        + q1y * (q1x * (s012 - s112) + q1w * (s022 - s122))
                        + qperpw * qperpy * s122
                        + qperpx * qperpz * s122
                        + q1z * (-(q1w * s012) + q1x * s022 + q1w * s112 - q1x * s122))
                        * theta,
                );

                c1[1] = derivative_term(
                    -t0y + t1y,
                    -(qperpx * qperpy * s000) - qperpw * qperpz * s000 - s010
                        + q1z * q1z * s010
                        + qperpx * qperpx * s010
                        + qperpz * qperpz * s010
                        - q1y * q1z * s020
                        + qperpw * qperpx * s020
                        - qperpy * qperpz * s020
                        + qperpx * qperpy * s100
                        + qperpw * qperpz * s100
                        + q1w * q1z * (-s000 + s100)
                        + q1x * q1x * (s010 - s110)
                        + s110
                        - q1z * q1z * s110
                        - qperpx * qperpx * s110
                        - qperpz * qperpz * s110
                        + q1x * (q1y * (-s000 + s100) + q1w * (s020 - s120))
                        + q1y * q1z * s120
                        - qperpw * qperpx * s120
                        + qperpy * qperpz * s120,
                    -(qperpx * qperpy * s001) - qperpw * qperpz * s001 - s011
                        + q1z * q1z * s011
                        + qperpx * qperpx * s011
                        + qperpz * qperpz * s011
                        - q1y * q1z * s021
                        + qperpw * qperpx * s021
                        - qperpy * qperpz * s021
                        + qperpx * qperpy * s101
                        + qperpw * qperpz * s101
                        + q1w * q1z * (-s001 + s101)
                        + q1x * q1x * (s011 - s111)
                        + s111
                        - q1z * q1z * s111
                        - qperpx * qperpx * s111
                        - qperpz * qperpz * s111
                        + q1x * (q1y * (-s001 + s101) + q1w * (s021 - s121))
                        + q1y * q1z * s121
                        - qperpw * qperpx * s121
                        + qperpy * qperpz * s121,
                    -(qperpx * qperpy * s002) - qperpw * qperpz * s002 - s012
                        + q1z * q1z * s012
                        + qperpx * qperpx * s012
                        + qperpz * qperpz * s012
                        - q1y * q1z * s022
                        + qperpw * qperpx * s022
                        - qperpy * qperpz * s022
                        + qperpx * qperpy * s102
                        + qperpw * qperpz * s102
                        + q1w * q1z * (-s002 + s102)
                        + q1x * q1x * (s012 - s112)
                        + s112
                        - q1z * q1z * s112
                        - qperpx * qperpx * s112
                        - qperpz * qperpz * s112
                        + q1x * (q1y * (-s002 + s102) + q1w * (s022 - s122))
                        + q1y * q1z * s122
                        - qperpw * qperpx * s122
                        + qperpy * qperpz * s122,
                );

                c2[1] = derivative_term(
                    0.0,
                    qperpx * qperpy * s000 + qperpw * qperpz * s000 + q1z * q1z * s010
                        - qperpx * qperpx * s010
                        - qperpz * qperpz * s010
                        - q1y * q1z * s020
                        - qperpw * qperpx * s020
                        + qperpy * qperpz * s020
                        - qperpx * qperpy * s100
                        - qperpw * qperpz * s100
                        + q1x * q1x * (s010 - s110)
                        - q1z * q1z * s110
                        + qperpx * qperpx * s110
                        + qperpz * qperpz * s110
                        + q1y * q1z * s120
                        + qperpw * qperpx * s120
                        - qperpy * qperpz * s120
                        + 2.0 * q1z * qperpw * s000 * theta
                        + 2.0 * q1y * qperpx * s000 * theta
                        - 4.0 * q1z * qperpz * s010 * theta
                        + 2.0 * q1z * qperpy * s020 * theta
                        + 2.0 * q1y * qperpz * s020 * theta
                        + q1x
                            * (q1w * s020 + q1y * (-s000 + s100) - q1w * s120
                                + 2.0 * qperpy * s000 * theta
                                - 4.0 * qperpx * s010 * theta
                                - 2.0 * qperpw * s020 * theta)
                        + q1w
                            * (-(q1z * s000) + q1z * s100 + 2.0 * qperpz * s000 * theta
                                - 2.0 * qperpx * s020 * theta),
                    qperpx * qperpy * s001 + qperpw * qperpz * s001 + q1z * q1z * s011
                        - qperpx * qperpx * s011
                        - qperpz * qperpz * s011
                        - q1y * q1z * s021
                        - qperpw * qperpx * s021
                        + qperpy * qperpz * s021
                        - qperpx * qperpy * s101
                        - qperpw * qperpz * s101
                        + q1x * q1x * (s011 - s111)
                        - q1z * q1z * s111
                        + qperpx * qperpx * s111
                        + qperpz * qperpz * s111
                        + q1y * q1z * s121
                        + qperpw * qperpx * s121
                        - qperpy * qperpz * s121
                        + 2.0 * q1z * qperpw * s001 * theta
                        + 2.0 * q1y * qperpx * s001 * theta
                        - 4.0 * q1z * qperpz * s011 * theta
                        + 2.0 * q1z * qperpy * s021 * theta
                        + 2.0 * q1y * qperpz * s021 * theta
                        + q1x
                            * (q1w * s021 + q1y * (-s001 + s101) - q1w * s121
                                + 2.0 * qperpy * s001 * theta
                                - 4.0 * qperpx * s011 * theta
                                - 2.0 * qperpw * s021 * theta)
                        + q1w
                            * (-(q1z * s001) + q1z * s101 + 2.0 * qperpz * s001 * theta
                                - 2.0 * qperpx * s021 * theta),
                    qperpx * qperpy * s002 + qperpw * qperpz * s002 + q1z * q1z * s012
                        - qperpx * qperpx * s012
                        - qperpz * qperpz * s012
                        - q1y * q1z * s022
                        - qperpw * qperpx * s022
                        + qperpy * qperpz * s022
                        - qperpx * qperpy * s102
                        - qperpw * qperpz * s102
                        + q1x * q1x * (s012 - s112)
                        - q1z * q1z * s112
                        + qperpx * qperpx * s112
                        + qperpz * qperpz * s112
                        + q1y * q1z * s122
                        + qperpw * qperpx * s122
                        - qperpy * qperpz * s122
                        + 2.0 * q1z * qperpw * s002 * theta
                        + 2.0 * q1y * qperpx * s002 * theta
                        - 4.0 * q1z * qperpz * s012 * theta
                        + 2.0 * q1z * qperpy * s022 * theta
                        + 2.0 * q1y * qperpz * s022 * theta
                        + q1x
                            * (q1w * s022 + q1y * (-s002 + s102) - q1w * s122
                                + 2.0 * qperpy * s002 * theta
                                - 4.0 * qperpx * s012 * theta
                                - 2.0 * qperpw * s022 * theta)
                        + q1w
                            * (-(q1z * s002) + q1z * s102 + 2.0 * qperpz * s002 * theta
                                - 2.0 * qperpx * s022 * theta),
                );

                c3[1] = derivative_term(
                    0.0,
                    2.0 * (-(q1x * qperpy * s000) - q1w * qperpz * s000
                        + 2.0 * q1x * qperpx * s010
                        + q1x * qperpw * s020
                        + q1w * qperpx * s020
                        + q1x * qperpy * s100
                        + q1w * qperpz * s100
                        - 2.0 * q1x * qperpx * s110
                        - q1x * qperpw * s120
                        - q1w * qperpx * s120
                        + q1z
                            * (2.0 * qperpz * s010 - qperpy * s020 + qperpw * (-s000 + s100)
                                - 2.0 * qperpz * s110
                                + qperpy * s120)
                        + q1y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 + qperpz * s120))
                        * theta,
                    2.0 * (-(q1x * qperpy * s001) - q1w * qperpz * s001
                        + 2.0 * q1x * qperpx * s011
                        + q1x * qperpw * s021
                        + q1w * qperpx * s021
                        + q1x * qperpy * s101
                        + q1w * qperpz * s101
                        - 2.0 * q1x * qperpx * s111
                        - q1x * qperpw * s121
                        - q1w * qperpx * s121
                        + q1z
                            * (2.0 * qperpz * s011 - qperpy * s021 + qperpw * (-s001 + s101)
                                - 2.0 * qperpz * s111
                                + qperpy * s121)
                        + q1y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 + qperpz * s121))
                        * theta,
                    2.0 * (-(q1x * qperpy * s002) - q1w * qperpz * s002
                        + 2.0 * q1x * qperpx * s012
                        + q1x * qperpw * s022
                        + q1w * qperpx * s022
                        + q1x * qperpy * s102
                        + q1w * qperpz * s102
                        - 2.0 * q1x * qperpx * s112
                        - q1x * qperpw * s122
                        - q1w * qperpx * s122
                        + q1z
                            * (2.0 * qperpz * s012 - qperpy * s022 + qperpw * (-s002 + s102)
                                - 2.0 * qperpz * s112
                                + qperpy * s122)
                        + q1y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 + qperpz * s122))
                        * theta,
                );

                c4[1] = derivative_term(
                    0.0,
                    -(q1x * qperpy * s000) - q1w * qperpz * s000
                        + 2.0 * q1x * qperpx * s010
                        + q1x * qperpw * s020
                        + q1w * qperpx * s020
                        + q1x * qperpy * s100
                        + q1w * qperpz * s100
                        - 2.0 * q1x * qperpx * s110
                        - q1x * qperpw * s120
                        - q1w * qperpx * s120
                        + 2.0 * qperpx * qperpy * s000 * theta
                        + 2.0 * qperpw * qperpz * s000 * theta
                        + 2.0 * q1x * q1x * s010 * theta
                        + 2.0 * q1z * q1z * s010 * theta
                        - 2.0 * qperpx * qperpx * s010 * theta
                        - 2.0 * qperpz * qperpz * s010 * theta
                        + 2.0 * q1w * q1x * s020 * theta
                        - 2.0 * qperpw * qperpx * s020 * theta
                        + 2.0 * qperpy * qperpz * s020 * theta
                        + q1y
                            * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 + qperpz * s120
                                - 2.0 * q1x * s000 * theta)
                        + q1z
                            * (2.0 * qperpz * s010 - qperpy * s020 + qperpw * (-s000 + s100)
                                - 2.0 * qperpz * s110
                                + qperpy * s120
                                - 2.0 * q1w * s000 * theta
                                - 2.0 * q1y * s020 * theta),
                    -(q1x * qperpy * s001) - q1w * qperpz * s001
                        + 2.0 * q1x * qperpx * s011
                        + q1x * qperpw * s021
                        + q1w * qperpx * s021
                        + q1x * qperpy * s101
                        + q1w * qperpz * s101
                        - 2.0 * q1x * qperpx * s111
                        - q1x * qperpw * s121
                        - q1w * qperpx * s121
                        + 2.0 * qperpx * qperpy * s001 * theta
                        + 2.0 * qperpw * qperpz * s001 * theta
                        + 2.0 * q1x * q1x * s011 * theta
                        + 2.0 * q1z * q1z * s011 * theta
                        - 2.0 * qperpx * qperpx * s011 * theta
                        - 2.0 * qperpz * qperpz * s011 * theta
                        + 2.0 * q1w * q1x * s021 * theta
                        - 2.0 * qperpw * qperpx * s021 * theta
                        + 2.0 * qperpy * qperpz * s021 * theta
                        + q1y
                            * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 + qperpz * s121
                                - 2.0 * q1x * s001 * theta)
                        + q1z
                            * (2.0 * qperpz * s011 - qperpy * s021 + qperpw * (-s001 + s101)
                                - 2.0 * qperpz * s111
                                + qperpy * s121
                                - 2.0 * q1w * s001 * theta
                                - 2.0 * q1y * s021 * theta),
                    -(q1x * qperpy * s002) - q1w * qperpz * s002
                        + 2.0 * q1x * qperpx * s012
                        + q1x * qperpw * s022
                        + q1w * qperpx * s022
                        + q1x * qperpy * s102
                        + q1w * qperpz * s102
                        - 2.0 * q1x * qperpx * s112
                        - q1x * qperpw * s122
                        - q1w * qperpx * s122
                        + 2.0 * qperpx * qperpy * s002 * theta
                        + 2.0 * qperpw * qperpz * s002 * theta
                        + 2.0 * q1x * q1x * s012 * theta
                        + 2.0 * q1z * q1z * s012 * theta
                        - 2.0 * qperpx * qperpx * s012 * theta
                        - 2.0 * qperpz * qperpz * s012 * theta
                        + 2.0 * q1w * q1x * s022 * theta
                        - 2.0 * qperpw * qperpx * s022 * theta
                        + 2.0 * qperpy * qperpz * s022 * theta
                        + q1y
                            * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 + qperpz * s122
                                - 2.0 * q1x * s002 * theta)
                        + q1z
                            * (2.0 * qperpz * s012 - qperpy * s022 + qperpw * (-s002 + s102)
                                - 2.0 * qperpz * s112
                                + qperpy * s122
                                - 2.0 * q1w * s002 * theta
                                - 2.0 * q1y * s022 * theta),
                );

                c5[1] = derivative_term(
                    0.0,
                    -2.0 * (qperpx * qperpy * s000 + qperpw * qperpz * s000 + q1z * q1z * s010
                        - qperpx * qperpx * s010
                        - qperpz * qperpz * s010
                        - q1y * q1z * s020
                        - qperpw * qperpx * s020
                        + qperpy * qperpz * s020
                        - qperpx * qperpy * s100
                        - qperpw * qperpz * s100
                        + q1w * q1z * (-s000 + s100)
                        + q1x * q1x * (s010 - s110)
                        - q1z * q1z * s110
                        + qperpx * qperpx * s110
                        + qperpz * qperpz * s110
                        + q1x * (q1y * (-s000 + s100) + q1w * (s020 - s120))
                        + q1y * q1z * s120
                        + qperpw * qperpx * s120
                        - qperpy * qperpz * s120)
                        * theta,
                    -2.0 * (qperpx * qperpy * s001 + qperpw * qperpz * s001 + q1z * q1z * s011
                        - qperpx * qperpx * s011
                        - qperpz * qperpz * s011
                        - q1y * q1z * s021
                        - qperpw * qperpx * s021
                        + qperpy * qperpz * s021
                        - qperpx * qperpy * s101
                        - qperpw * qperpz * s101
                        + q1w * q1z * (-s001 + s101)
                        + q1x * q1x * (s011 - s111)
                        - q1z * q1z * s111
                        + qperpx * qperpx * s111
                        + qperpz * qperpz * s111
                        + q1x * (q1y * (-s001 + s101) + q1w * (s021 - s121))
                        + q1y * q1z * s121
                        + qperpw * qperpx * s121
                        - qperpy * qperpz * s121)
                        * theta,
                    -2.0 * (qperpx * qperpy * s002 + qperpw * qperpz * s002 + q1z * q1z * s012
                        - qperpx * qperpx * s012
                        - qperpz * qperpz * s012
                        - q1y * q1z * s022
                        - qperpw * qperpx * s022
                        + qperpy * qperpz * s022
                        - qperpx * qperpy * s102
                        - qperpw * qperpz * s102
                        + q1w * q1z * (-s002 + s102)
                        + q1x * q1x * (s012 - s112)
                        - q1z * q1z * s112
                        + qperpx * qperpx * s112
                        + qperpz * qperpz * s112
                        + q1x * (q1y * (-s002 + s102) + q1w * (s022 - s122))
                        + q1y * q1z * s122
                        + qperpw * qperpx * s122
                        - qperpy * qperpz * s122)
                        * theta,
                );

                c1[2] = derivative_term(
                    -t0z + t1z,
                    qperpw * qperpy * s000
                        - qperpx * qperpz * s000
                        - q1y * q1z * s010
                        - qperpw * qperpx * s010
                        - qperpy * qperpz * s010
                        - s020
                        + q1y * q1y * s020
                        + qperpx * qperpx * s020
                        + qperpy * qperpy * s020
                        - qperpw * qperpy * s100
                        + qperpx * qperpz * s100
                        + q1x * q1z * (-s000 + s100)
                        + q1y * q1z * s110
                        + qperpw * qperpx * s110
                        + qperpy * qperpz * s110
                        + q1w * (q1y * (s000 - s100) + q1x * (-s010 + s110))
                        + q1x * q1x * (s020 - s120)
                        + s120
                        - q1y * q1y * s120
                        - qperpx * qperpx * s120
                        - qperpy * qperpy * s120,
                    qperpw * qperpy * s001
                        - qperpx * qperpz * s001
                        - q1y * q1z * s011
                        - qperpw * qperpx * s011
                        - qperpy * qperpz * s011
                        - s021
                        + q1y * q1y * s021
                        + qperpx * qperpx * s021
                        + qperpy * qperpy * s021
                        - qperpw * qperpy * s101
                        + qperpx * qperpz * s101
                        + q1x * q1z * (-s001 + s101)
                        + q1y * q1z * s111
                        + qperpw * qperpx * s111
                        + qperpy * qperpz * s111
                        + q1w * (q1y * (s001 - s101) + q1x * (-s011 + s111))
                        + q1x * q1x * (s021 - s121)
                        + s121
                        - q1y * q1y * s121
                        - qperpx * qperpx * s121
                        - qperpy * qperpy * s121,
                    qperpw * qperpy * s002
                        - qperpx * qperpz * s002
                        - q1y * q1z * s012
                        - qperpw * qperpx * s012
                        - qperpy * qperpz * s012
                        - s022
                        + q1y * q1y * s022
                        + qperpx * qperpx * s022
                        + qperpy * qperpy * s022
                        - qperpw * qperpy * s102
                        + qperpx * qperpz * s102
                        + q1x * q1z * (-s002 + s102)
                        + q1y * q1z * s112
                        + qperpw * qperpx * s112
                        + qperpy * qperpz * s112
                        + q1w * (q1y * (s002 - s102) + q1x * (-s012 + s112))
                        + q1x * q1x * (s022 - s122)
                        + s122
                        - q1y * q1y * s122
                        - qperpx * qperpx * s122
                        - qperpy * qperpy * s122,
                );

                c2[2] = derivative_term(
                    0.0,
                    q1w * q1y * s000 - q1x * q1z * s000 - qperpw * qperpy * s000
                        + qperpx * qperpz * s000
                        - q1w * q1x * s010
                        - q1y * q1z * s010
                        + qperpw * qperpx * s010
                        + qperpy * qperpz * s010
                        + q1x * q1x * s020
                        + q1y * q1y * s020
                        - qperpx * qperpx * s020
                        - qperpy * qperpy * s020
                        - q1w * q1y * s100
                        + q1x * q1z * s100
                        + qperpw * qperpy * s100
                        - qperpx * qperpz * s100
                        + q1w * q1x * s110
                        + q1y * q1z * s110
                        - qperpw * qperpx * s110
                        - qperpy * qperpz * s110
                        - q1x * q1x * s120
                        - q1y * q1y * s120
                        + qperpx * qperpx * s120
                        + qperpy * qperpy * s120
                        - 2.0 * q1y * qperpw * s000 * theta
                        + 2.0 * q1z * qperpx * s000 * theta
                        - 2.0 * q1w * qperpy * s000 * theta
                        + 2.0 * q1x * qperpz * s000 * theta
                        + 2.0 * q1x * qperpw * s010 * theta
                        + 2.0 * q1w * qperpx * s010 * theta
                        + 2.0 * q1z * qperpy * s010 * theta
                        + 2.0 * q1y * qperpz * s010 * theta
                        - 4.0 * q1x * qperpx * s020 * theta
                        - 4.0 * q1y * qperpy * s020 * theta,
                    q1w * q1y * s001 - q1x * q1z * s001 - qperpw * qperpy * s001
                        + qperpx * qperpz * s001
                        - q1w * q1x * s011
                        - q1y * q1z * s011
                        + qperpw * qperpx * s011
                        + qperpy * qperpz * s011
                        + q1x * q1x * s021
                        + q1y * q1y * s021
                        - qperpx * qperpx * s021
                        - qperpy * qperpy * s021
                        - q1w * q1y * s101
                        + q1x * q1z * s101
                        + qperpw * qperpy * s101
                        - qperpx * qperpz * s101
                        + q1w * q1x * s111
                        + q1y * q1z * s111
                        - qperpw * qperpx * s111
                        - qperpy * qperpz * s111
                        - q1x * q1x * s121
                        - q1y * q1y * s121
                        + qperpx * qperpx * s121
                        + qperpy * qperpy * s121
                        - 2.0 * q1y * qperpw * s001 * theta
                        + 2.0 * q1z * qperpx * s001 * theta
                        - 2.0 * q1w * qperpy * s001 * theta
                        + 2.0 * q1x * qperpz * s001 * theta
                        + 2.0 * q1x * qperpw * s011 * theta
                        + 2.0 * q1w * qperpx * s011 * theta
                        + 2.0 * q1z * qperpy * s011 * theta
                        + 2.0 * q1y * qperpz * s011 * theta
                        - 4.0 * q1x * qperpx * s021 * theta
                        - 4.0 * q1y * qperpy * s021 * theta,
                    q1w * q1y * s002 - q1x * q1z * s002 - qperpw * qperpy * s002
                        + qperpx * qperpz * s002
                        - q1w * q1x * s012
                        - q1y * q1z * s012
                        + qperpw * qperpx * s012
                        + qperpy * qperpz * s012
                        + q1x * q1x * s022
                        + q1y * q1y * s022
                        - qperpx * qperpx * s022
                        - qperpy * qperpy * s022
                        - q1w * q1y * s102
                        + q1x * q1z * s102
                        + qperpw * qperpy * s102
                        - qperpx * qperpz * s102
                        + q1w * q1x * s112
                        + q1y * q1z * s112
                        - qperpw * qperpx * s112
                        - qperpy * qperpz * s112
                        - q1x * q1x * s122
                        - q1y * q1y * s122
                        + qperpx * qperpx * s122
                        + qperpy * qperpy * s122
                        - 2.0 * q1y * qperpw * s002 * theta
                        + 2.0 * q1z * qperpx * s002 * theta
                        - 2.0 * q1w * qperpy * s002 * theta
                        + 2.0 * q1x * qperpz * s002 * theta
                        + 2.0 * q1x * qperpw * s012 * theta
                        + 2.0 * q1w * qperpx * s012 * theta
                        + 2.0 * q1z * qperpy * s012 * theta
                        + 2.0 * q1y * qperpz * s012 * theta
                        - 4.0 * q1x * qperpx * s022 * theta
                        - 4.0 * q1y * qperpy * s022 * theta,
                );

                c3[2] = derivative_term(
                    0.0,
                    -2.0 * (-(q1w * qperpy * s000)
                        + q1x * qperpz * s000
                        + q1x * qperpw * s010
                        + q1w * qperpx * s010
                        - 2.0 * q1x * qperpx * s020
                        + q1w * qperpy * s100
                        - q1x * qperpz * s100
                        - q1x * qperpw * s110
                        - q1w * qperpx * s110
                        + q1z * (qperpx * s000 + qperpy * s010 - qperpx * s100 - qperpy * s110)
                        + 2.0 * q1x * qperpx * s120
                        + q1y
                            * (qperpz * s010 - 2.0 * qperpy * s020 + qperpw * (-s000 + s100)
                                - qperpz * s110
                                + 2.0 * qperpy * s120))
                        * theta,
                    -2.0 * (-(q1w * qperpy * s001)
                        + q1x * qperpz * s001
                        + q1x * qperpw * s011
                        + q1w * qperpx * s011
                        - 2.0 * q1x * qperpx * s021
                        + q1w * qperpy * s101
                        - q1x * qperpz * s101
                        - q1x * qperpw * s111
                        - q1w * qperpx * s111
                        + q1z * (qperpx * s001 + qperpy * s011 - qperpx * s101 - qperpy * s111)
                        + 2.0 * q1x * qperpx * s121
                        + q1y
                            * (qperpz * s011 - 2.0 * qperpy * s021 + qperpw * (-s001 + s101)
                                - qperpz * s111
                                + 2.0 * qperpy * s121))
                        * theta,
                    -2.0 * (-(q1w * qperpy * s002)
                        + q1x * qperpz * s002
                        + q1x * qperpw * s012
                        + q1w * qperpx * s012
                        - 2.0 * q1x * qperpx * s022
                        + q1w * qperpy * s102
                        - q1x * qperpz * s102
                        - q1x * qperpw * s112
                        - q1w * qperpx * s112
                        + q1z * (qperpx * s002 + qperpy * s012 - qperpx * s102 - qperpy * s112)
                        + 2.0 * q1x * qperpx * s122
                        + q1y
                            * (qperpz * s012 - 2.0 * qperpy * s022 + qperpw * (-s002 + s102)
                                - qperpz * s112
                                + 2.0 * qperpy * s122))
                        * theta,
                );

                c4[2] = derivative_term(
                    0.0,
                    q1w * qperpy * s000
                        - q1x * qperpz * s000
                        - q1x * qperpw * s010
                        - q1w * qperpx * s010
                        + 2.0 * q1x * qperpx * s020
                        - q1w * qperpy * s100
                        + q1x * qperpz * s100
                        + q1x * qperpw * s110
                        + q1w * qperpx * s110
                        - 2.0 * q1x * qperpx * s120
                        - 2.0 * qperpw * qperpy * s000 * theta
                        + 2.0 * qperpx * qperpz * s000 * theta
                        - 2.0 * q1w * q1x * s010 * theta
                        + 2.0 * qperpw * qperpx * s010 * theta
                        + 2.0 * qperpy * qperpz * s010 * theta
                        + 2.0 * q1x * q1x * s020 * theta
                        + 2.0 * q1y * q1y * s020 * theta
                        - 2.0 * qperpx * qperpx * s020 * theta
                        - 2.0 * qperpy * qperpy * s020 * theta
                        + q1z
                            * (-(qperpx * s000) - qperpy * s010 + qperpx * s100 + qperpy * s110
                                - 2.0 * q1x * s000 * theta)
                        + q1y
                            * (-(qperpz * s010)
                                + 2.0 * qperpy * s020
                                + qperpw * (s000 - s100)
                                + qperpz * s110
                                - 2.0 * qperpy * s120
                                + 2.0 * q1w * s000 * theta
                                - 2.0 * q1z * s010 * theta),
                    q1w * qperpy * s001
                        - q1x * qperpz * s001
                        - q1x * qperpw * s011
                        - q1w * qperpx * s011
                        + 2.0 * q1x * qperpx * s021
                        - q1w * qperpy * s101
                        + q1x * qperpz * s101
                        + q1x * qperpw * s111
                        + q1w * qperpx * s111
                        - 2.0 * q1x * qperpx * s121
                        - 2.0 * qperpw * qperpy * s001 * theta
                        + 2.0 * qperpx * qperpz * s001 * theta
                        - 2.0 * q1w * q1x * s011 * theta
                        + 2.0 * qperpw * qperpx * s011 * theta
                        + 2.0 * qperpy * qperpz * s011 * theta
                        + 2.0 * q1x * q1x * s021 * theta
                        + 2.0 * q1y * q1y * s021 * theta
                        - 2.0 * qperpx * qperpx * s021 * theta
                        - 2.0 * qperpy * qperpy * s021 * theta
                        + q1z
                            * (-(qperpx * s001) - qperpy * s011 + qperpx * s101 + qperpy * s111
                                - 2.0 * q1x * s001 * theta)
                        + q1y
                            * (-(qperpz * s011)
                                + 2.0 * qperpy * s021
                                + qperpw * (s001 - s101)
                                + qperpz * s111
                                - 2.0 * qperpy * s121
                                + 2.0 * q1w * s001 * theta
                                - 2.0 * q1z * s011 * theta),
                    q1w * qperpy * s002
                        - q1x * qperpz * s002
                        - q1x * qperpw * s012
                        - q1w * qperpx * s012
                        + 2.0 * q1x * qperpx * s022
                        - q1w * qperpy * s102
                        + q1x * qperpz * s102
                        + q1x * qperpw * s112
                        + q1w * qperpx * s112
                        - 2.0 * q1x * qperpx * s122
                        - 2.0 * qperpw * qperpy * s002 * theta
                        + 2.0 * qperpx * qperpz * s002 * theta
                        - 2.0 * q1w * q1x * s012 * theta
                        + 2.0 * qperpw * qperpx * s012 * theta
                        + 2.0 * qperpy * qperpz * s012 * theta
                        + 2.0 * q1x * q1x * s022 * theta
                        + 2.0 * q1y * q1y * s022 * theta
                        - 2.0 * qperpx * qperpx * s022 * theta
                        - 2.0 * qperpy * qperpy * s022 * theta
                        + q1z
                            * (-(qperpx * s002) - qperpy * s012 + qperpx * s102 + qperpy * s112
                                - 2.0 * q1x * s002 * theta)
                        + q1y
                            * (-(qperpz * s012)
                                + 2.0 * qperpy * s022
                                + qperpw * (s002 - s102)
                                + qperpz * s112
                                - 2.0 * qperpy * s122
                                + 2.0 * q1w * s002 * theta
                                - 2.0 * q1z * s012 * theta),
                );

                c5[2] = derivative_term(
                    0.0,
                    2.0 * (qperpw * qperpy * s000 - qperpx * qperpz * s000 + q1y * q1z * s010
                        - qperpw * qperpx * s010
                        - qperpy * qperpz * s010
                        - q1y * q1y * s020
                        + qperpx * qperpx * s020
                        + qperpy * qperpy * s020
                        + q1x * q1z * (s000 - s100)
                        - qperpw * qperpy * s100
                        + qperpx * qperpz * s100
                        + q1w * (q1y * (-s000 + s100) + q1x * (s010 - s110))
                        - q1y * q1z * s110
                        + qperpw * qperpx * s110
                        + qperpy * qperpz * s110
                        + q1y * q1y * s120
                        - qperpx * qperpx * s120
                        - qperpy * qperpy * s120
                        + q1x * q1x * (-s020 + s120))
                        * theta,
                    2.0 * (qperpw * qperpy * s001 - qperpx * qperpz * s001 + q1y * q1z * s011
                        - qperpw * qperpx * s011
                        - qperpy * qperpz * s011
                        - q1y * q1y * s021
                        + qperpx * qperpx * s021
                        + qperpy * qperpy * s021
                        + q1x * q1z * (s001 - s101)
                        - qperpw * qperpy * s101
                        + qperpx * qperpz * s101
                        + q1w * (q1y * (-s001 + s101) + q1x * (s011 - s111))
                        - q1y * q1z * s111
                        + qperpw * qperpx * s111
                        + qperpy * qperpz * s111
                        + q1y * q1y * s121
                        - qperpx * qperpx * s121
                        - qperpy * qperpy * s121
                        + q1x * q1x * (-s021 + s121))
                        * theta,
                    2.0 * (qperpw * qperpy * s002 - qperpx * qperpz * s002 + q1y * q1z * s012
                        - qperpw * qperpx * s012
                        - qperpy * qperpz * s012
                        - q1y * q1y * s022
                        + qperpx * qperpx * s022
                        + qperpy * qperpy * s022
                        + q1x * q1z * (s002 - s102)
                        - qperpw * qperpy * s102
                        + qperpx * qperpz * s102
                        + q1w * (q1y * (-s002 + s102) + q1x * (s012 - s112))
                        - q1y * q1z * s112
                        + qperpw * qperpx * s112
                        + qperpy * qperpz * s112
                        + q1y * q1y * s122
                        - qperpx * qperpx * s122
                        - qperpy * qperpy * s122
                        + q1x * q1x * (-s022 + s122))
                        * theta,
                );
            }
        }

        Self {
            start_transform,
            end_transform,
            start_time,
            end_time,
            actually_animated,
            t,
            r,
            s,
            has_rotation,
            c1,
            c2,
            c3,
            c4,
            c5,
        }
    }

    /// Returns the interpolated trnasformation matrix at the given time.
    ///
    /// * `time` - Time.
    pub fn interpolate(&self, time: Float) -> ArcTransform {
        // Handle boundary conditions for matrix interpolation.
        if !self.actually_animated || time <= self.start_time {
            return Arc::clone(&self.start_transform);
        }
        if time >= self.end_time {
            return Arc::clone(&self.end_transform);
        }

        // Map time to range [0, 1].
        let dt = (time - self.start_time) / (self.end_time - self.start_time);

        // Interpolate translation at dt.
        let trans = (1.0 - dt) * self.t[0] + dt * self.t[1];

        // Interpolate rotation at dt.
        let rotate = self.r[0].slerp(dt, self.r[1]);

        // Interpolate scale at dt.
        let mut scale = Matrix4x4::default();
        for i in 0..3 {
            for j in 0..3 {
                scale.m[i][j] = lerp(dt, self.s[0].m[i][j], self.s[1].m[i][j]);
            }
        }

        // Compute interpolated matrix as product of interpolated components
        Arc::new(Transform::translate(&trans) * Transform::from(rotate) * Transform::from(scale))
    }

    /// Applies animated transformation to a given ray.
    ///
    /// * `r` - The ray.
    pub fn transform_ray(&self, r: &Ray) -> Ray {
        let t = if !self.actually_animated || r.time <= self.start_time {
            Arc::clone(&self.start_transform)
        } else if r.time >= self.end_time {
            Arc::clone(&self.end_transform)
        } else {
            self.interpolate(r.time)
        };
        t.transform_ray(r)
    }

    /// Applies animated transformation to a given point.
    ///
    /// * `time` - The time.
    /// * `p`    - The point.
    pub fn transform_point(&self, time: Float, p: &Point3f) -> Point3f {
        let t = if !self.actually_animated || time <= self.start_time {
            Arc::clone(&self.start_transform)
        } else if time >= self.end_time {
            Arc::clone(&self.end_transform)
        } else {
            self.interpolate(time)
        };
        t.transform_point(p)
    }

    /// Applies animated transformation to a given vector.
    ///
    /// * `time` - The time.
    /// * `v`    - The vector.
    pub fn transform_vector(&self, time: Float, v: &Vector3f) -> Vector3f {
        let t = if !self.actually_animated || time <= self.start_time {
            Arc::clone(&self.start_transform)
        } else if time >= self.end_time {
            Arc::clone(&self.end_transform)
        } else {
            self.interpolate(time)
        };
        t.transform_vector(v)
    }

    /// Compute a bounding box that encompasses all of its motion over the
    /// animation time period.
    ///
    /// * `b` - The bounding box.
    pub fn motion_bounds(&self, b: &Bounds3f) -> Bounds3f {
        if !self.actually_animated {
            return self.start_transform.transform_bounds(b);
        }

        if !self.has_rotation {
            return self
                .start_transform
                .transform_bounds(b)
                .union(&self.end_transform.transform_bounds(b));
        }

        (0..8u8).fold(Bounds3::<Float>::empty(), |bounds, corner| {
            bounds.union(&self.bound_point_motion(&b.corner(corner)))
        })
    }

    /// Use the coefficients `ki`, based on the keyframe matrices to compute
    /// a robust bound of the motion of point `p`.
    ///
    /// * `p` - The point.
    pub fn bound_point_motion(&self, p: &Point3f) -> Bounds3f {
        if !self.actually_animated {
            return Bounds3::from(self.start_transform.transform_point(p));
        }

        let mut bounds = Bounds3::new(
            self.start_transform.transform_point(p),
            self.end_transform.transform_point(p),
        );

        let cos_theta = self.r[0].dot(&self.r[1]);
        let theta = clamp(cos_theta, -1.0, 1.0).acos();
        for c in 0..3 {
            // Find any motion derivative zeros for the component c.
            let mut zeros = [0.0; 8];
            let mut n_zeros = 0_usize;
            Interval::new(0.0, 1.0).find_zeros(
                self.c1[c].eval(p),
                self.c2[c].eval(p),
                self.c3[c].eval(p),
                self.c4[c].eval(p),
                self.c5[c].eval(p),
                theta,
                &mut zeros,
                &mut n_zeros,
                8, // zeros.len()
            );

            // Expand bounding box for any motion derivative zeros found
            bounds = (0..n_zeros).fold(bounds, |b, i| {
                let pz = self.transform_point(lerp(zeros[i], self.start_time, self.end_time), p);
                b.union(&pz)
            });
        }

        bounds
    }
}

/// Decomposes a transformation matrix into its translation, rotation and
/// scaling components.
///
/// The transformation decomposition is: M = TRS.
///
/// * `m`      - The matrix to decompose.
/// * `t`      - The translation component.
/// * `r_quat` - The rotation component.
/// * `s`      - The scaling component.
fn decompose(m: &Matrix4x4, t: &mut Vector3f, r_quat: &mut Quaternion, s: &mut Matrix4x4) {
    // Extract translation T from transformation matrix
    *t = Vector3f::new(m[0][3], m[1][3], m[2][3]);

    // Compute new transformation matrix M without translation
    let mut m1 = *m;
    for i in 0..3 {
        m1.m[i][3] = 0.0;
        m1.m[3][i] = 0.0;
    }
    m1.m[3][3] = 1.0;

    // Extract rotation R from transformation matrix
    let mut r = m1;
    let mut norm = 0.0;
    let mut count = 0;
    loop {
        // Compute the next matrix R_next in series
        let mut r_next = Matrix4x4::default();
        let r_it = r.transpose().inverse();
        for i in 0..4 {
            for j in 0..4 {
                r_next.m[i][j] = 0.5 * (r[i][j] + r_it[i][j]);
            }
        }

        // Compute norm of difference between R and R_next
        for i in 0..3 {
            let n = abs(r[i][0] - r_next[i][0])
                + abs(r[i][1] - r_next[i][1])
                + abs(r[i][2] - r_next[i][2]);
            norm = max(norm, n);
        }

        r = r_next;

        count += 1;
        if count >= 100 && norm <= 0.0001 {
            break;
        }
    }
    *r_quat = Quaternion::from(Transform::from(r));

    // Compute scale S using rotation and original matrix
    *s = r.inverse() * *m;
}

/// DerivativeTerm encapsulates the coefficients `ki` to bound the motion of a
/// point `p`.
#[derive(Copy, Clone, Debug, Default)]
struct DerivativeTerm {
    pub kc: Float,
    pub kx: Float,
    pub ky: Float,
    pub kz: Float,
}

/// Returns a new derivative term with given coefficients `ki`.
fn derivative_term(kc: Float, kx: Float, ky: Float, kz: Float) -> DerivativeTerm {
    DerivativeTerm { kc, kx, ky, kz }
}

impl DerivativeTerm {
    /// Return the coefficient `ci` to bound the motion of a point `p`.
    ///
    /// * `p` - The point.
    pub fn eval(&self, p: &Point3f) -> Float {
        self.kc + self.kx * p.x + self.ky * p.y + self.kz * p.z
    }
}
