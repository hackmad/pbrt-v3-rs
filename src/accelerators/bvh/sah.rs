//! Surface Area Heuristic Algorithm

#![allow(dead_code)]
use super::common::*;
use super::{ArcPrimitive, Axis, Bounds3f, Float, Union, N_BUCKETS};
use std::cmp::Ordering;
use std::sync::{Arc, Mutex};

/// Recursively build the BVH structure for either Middle, EqualCounts or SAH
/// algorithm.
///
/// * `primitives`        - The primitives in the node.
/// * `split_method`      - Middle|EqualCounts|SAH
/// * `max_prims_in_node` - Maximum number of primitives in the node.
/// * `primitive_info`    - Primitive information.
/// * `start`             - Starting index.
/// * `end`               - Ending index.
/// * `total_nodes`       - Used to return total number of nodes.
/// * `ordered_prims`     - Used to return a list of primitives ordered such that
///                         primitives in leaf nodes occupy contiguous ranges in
///                         the vector.
pub fn recursive_build(
    primitives: &Vec<ArcPrimitive>,
    split_method: SplitMethod,
    max_prims_in_node: u8,
    primitive_info: &mut Vec<BVHPrimitiveInfo>,
    start: usize,
    end: usize,
    total_nodes: &mut usize,
    ordered_prims: Arc<Mutex<Vec<ArcPrimitive>>>,
) -> Arc<BVHBuildNode> {
    // Compute bounds of all primitives in BVH node.
    let bounds = (start..end).fold(Bounds3f::default(), |b, i| {
        b.union(&primitive_info[i].bounds)
    });

    let mut dim = Axis::default(); // Will be set if we need to make interior node.

    let n_primitives = end - start;

    let interior_midpoint = if n_primitives == 1 {
        // Create leaf BVHBuildNode.
        None
    } else {
        // Compute bound of primitive centroids, choose split dimension dim.
        let centroid_bounds = (start..end).fold(Bounds3f::default(), |b, i| {
            b.union(&primitive_info[i].centroid)
        });
        dim = centroid_bounds.maximum_extent();

        // Partition primitives into two sets and build children.
        if centroid_bounds.p_max[dim] == centroid_bounds.p_min[dim] {
            // Create leaf BVHBuildNode.
            None
        } else {
            // Partition primitives based on split_method.
            match split_method {
                SplitMethod::Middle => Some(split_middle(
                    primitive_info,
                    start,
                    end,
                    dim,
                    &centroid_bounds,
                )),

                SplitMethod::EqualCounts => {
                    Some(split_equal_counts(primitive_info, start, end, dim))
                }

                SplitMethod::SAH => split_sah(
                    primitive_info,
                    start,
                    end,
                    n_primitives,
                    dim,
                    &centroid_bounds,
                    &bounds,
                    max_prims_in_node,
                ),

                _ => panic!("recursive_build(): Invalid split_method={:?}", split_method),
            }
        }
    };

    *total_nodes += 1;
    if let Some(mid) = interior_midpoint {
        // Create interior BVHBuildNode.
        let c0 = recursive_build(
            primitives,
            split_method,
            max_prims_in_node,
            primitive_info,
            start,
            mid,
            total_nodes,
            ordered_prims.clone(),
        );
        let c1 = recursive_build(
            primitives,
            split_method,
            max_prims_in_node,
            primitive_info,
            mid,
            end,
            total_nodes,
            ordered_prims.clone(),
        );
        create_bvh_interior_node(dim, c0, c1)
    } else {
        // Create leaf BVHBuildNode.
        let prims = ordered_prims.clone();
        let mut prims2 = prims.lock().expect("unable to lock ordered_prims");
        let first_prim_offset = prims2.len();
        for i in start..end {
            let prim_num = primitive_info[i].primitive_number;
            prims2.push(primitives[prim_num].clone());
        }
        create_bvh_leaf_node(first_prim_offset, n_primitives, bounds)
    }
}

/// Linear Bounding Volume Hierarchy using splitting planes that are
/// midpoint of each region of space.
///
/// * `primitive_info`  - Vector containing all primitive info.
/// * `start`           - Start index in primitive_info.
/// * `end`             - End index in primitive_info.
/// * `dim`             - Axis used to partition primitives.
/// * `centroid_bounds` - Bounding box of primtive centroids in primtive_info
///                       from start to end.
fn split_middle(
    primitive_info: &mut Vec<BVHPrimitiveInfo>,
    start: usize,
    end: usize,
    dim: Axis,
    centroid_bounds: &Bounds3f,
) -> usize {
    let pmid = (centroid_bounds.p_min[dim] + centroid_bounds.p_max[dim]) / 2.0;
    let split = primitive_info[start..end + 1]
        .iter_mut()
        .partition_in_place(|pi| pi.centroid[dim] < pmid);
    let mid = start + split;

    if mid != start && mid != end {
        mid
    } else {
        // Lots of prims with large overlapping bounding boxes, this may fail
        // to partition; in that case don't use EqualCounts.
        split_equal_counts(primitive_info, start, end, dim)
    }
}

/// Partition primitives into equally sized subsets such that the first half
/// of the primitives have smallest centroid coordinate values along the
/// chosen axis, and second have have the largest centroid coordinate values.
///
/// * `primitive_info`  - Vector containing all primitive info.
/// * `start`           - Start index in primitive_info.
/// * `end`             - End index in primitive_info.
/// * `dim`             - Axis used to partition primitives.
fn split_equal_counts(
    primitive_info: &mut Vec<BVHPrimitiveInfo>,
    start: usize,
    end: usize,
    dim: Axis,
) -> usize {
    let mid = (start + end) / 2;

    kth_element_by(primitive_info, start, mid, end, |&a, &b| {
        if a.centroid[dim] < b.centroid[dim] {
            Ordering::Less
        } else if a.centroid[dim] == b.centroid[dim] {
            Ordering::Equal
        } else {
            Ordering::Greater
        }
    });

    mid
}

/// Partition primitives using Surface Area Heuristic.
///
/// If the algorithm is able to partition primitives it will return the pivot
/// index (mid) for interior node creation; otherwise None is returned to
/// indicate leaf node creation.
///
/// * `primitive_info`    - Vector containing all primitive info.
/// * `start`             - Start index in primitive_info.
/// * `end`               - End index in primitive_info.
/// * `n_primitives`      - Number of primitives between start and end.
/// * `dim`               - Axis used to partition primitives.
/// * `centroid_bounds`   - Bounding box of primtive centroids in primtive_info
///                         from start to end.
/// * `bounds`            - Bound box of all primitives in BVH node.
/// * `max_prims_in_node` - Maximum primitives allowed in node.
fn split_sah(
    primitive_info: &mut Vec<BVHPrimitiveInfo>,
    start: usize,
    end: usize,
    n_primitives: usize,
    dim: Axis,
    centroid_bounds: &Bounds3f,
    bounds: &Bounds3f,
    max_prims_in_node: u8,
) -> Option<usize> {
    // Partition primitives using approximate SAH.
    if n_primitives <= 2 {
        // Partition primitives into equally-sized subsets.
        Some(split_equal_counts(primitive_info, start, end, dim))
    } else {
        // Allocate BucketInfo for SAH partition buckets.
        let mut buckets = [BucketInfo::default(); N_BUCKETS];

        // Initialize BucketInfo for SAH partition buckets.
        for i in start..end {
            let mut b = (N_BUCKETS as Float
                * centroid_bounds.offset(&primitive_info[i].centroid)[dim])
                as usize;
            if b == N_BUCKETS {
                b = N_BUCKETS - 1;
            }
            debug_assert!(b > 0);
            debug_assert!(b < N_BUCKETS);
            buckets[b].count += 1;
            buckets[b].bounds = buckets[b].bounds.union(&primitive_info[i].bounds);
        }

        // Compute costs for splitting after each bucket
        let mut cost = [0.0 as Float; N_BUCKETS - 1];
        for i in 0..N_BUCKETS - 1 {
            let (mut b0, mut b1) = (Bounds3f::default(), Bounds3f::default());
            let (mut count0, mut count1) = (0, 0);
            for j in 0..i + 1 {
                b0 = b0.union(&buckets[j].bounds);
                count0 += buckets[j].count;
            }
            for j in i + 1..N_BUCKETS {
                b1 = b1.union(&buckets[j].bounds);
                count1 += buckets[j].count;
            }
            cost[i] = 1.0
                + (count0 as Float * b0.surface_area() + count1 as Float * b1.surface_area())
                    / bounds.surface_area();
        }

        // Find bucket to split at that minimizes SAH metric.
        let mut min_cost = cost[0];
        let mut min_cost_split_bucket = 0;
        for i in 1..N_BUCKETS - 1 {
            if cost[i] < min_cost {
                min_cost = cost[i];
                min_cost_split_bucket = i;
            }
        }

        let leaf_cost = n_primitives as Float;
        if n_primitives > max_prims_in_node as usize || min_cost < leaf_cost {
            // Partition primitives at selected SAH bucket and return the
            // pivot point as mid.
            let split = primitive_info[start..end + 1]
                .iter_mut()
                .partition_in_place(|pi| {
                    let mut b =
                        (N_BUCKETS as Float * centroid_bounds.offset(&pi.centroid)[dim]) as usize;
                    if b == N_BUCKETS {
                        b = N_BUCKETS - 1;
                    }
                    debug_assert!(b > 0);
                    debug_assert!(b < N_BUCKETS);
                    b <= min_cost_split_bucket
                });
            let pmid = start + split;
            Some(pmid)
        } else {
            // No split occurred. Indicate creation of leaf node.
            None
        }
    }
}
