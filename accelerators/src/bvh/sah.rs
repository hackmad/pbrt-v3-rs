//! Surface Area Heuristic Algorithm

use super::common::*;
use core::geometry::*;
use core::pbrt::*;
use core::primitive::*;
use order_stat::kth_by;
use shared_arena::{ArenaArc, SharedArena};
use std::cmp::Ordering;
use std::sync::{Arc, Mutex};

const N_BUCKETS: usize = 12;

/// Recursively build the BVH structure for either Middle, EqualCounts or SAH algorithm.
///
/// * `arena`             - Shared arena for memory allocations.
/// * `primitives`        - The primitives in the node.
/// * `split_method`      - Middle|EqualCounts|SAH
/// * `max_prims_in_node` - Maximum number of primitives in the node.
/// * `primitive_info`    - Primitive information.
/// * `start`             - Starting index. For first call it should be 0.
/// * `end`               - Ending index + 1. For first call it should be number of primitives.
/// * `total_nodes`       - Used to return total number of nodes.
/// * `ordered_prims`     - Used to return a list of primitives ordered such that primitives in leaf nodes occupy
///                         contiguous ranges in the vector.
pub fn build(
    arena: &SharedArena<BVHBuildNode>,
    primitives: &[ArcPrimitive],
    split_method: SplitMethod,
    max_prims_in_node: u8,
    primitive_info: &mut [BVHPrimitiveInfo],
    start: usize,
    end: usize,
    total_nodes: &mut usize,
    ordered_prims: Arc<Mutex<Vec<ArcPrimitive>>>,
) -> ArenaArc<BVHBuildNode> {
    assert_ne!(start, end);

    *total_nodes += 1;

    // Compute bounds of all primitives in BVH node.
    let mut bounds = Bounds3f::EMPTY;
    for i in start..end {
        bounds = bounds.union(&primitive_info[i].bounds);
    }

    let n_primitives = end - start;

    if n_primitives == 1 {
        // Create leaf BVHBuildNode.
        new_leaf_node(arena, primitives, primitive_info, start, end, ordered_prims, bounds)
    } else {
        // Compute bound of primitive centroids, choose split dimension dim.
        let mut centroid_bounds = Bounds3f::EMPTY;
        for i in start..end {
            centroid_bounds = centroid_bounds.union(&primitive_info[i].centroid);
        }
        let dim = centroid_bounds.maximum_extent();

        // Partition primitives into two sets and build children.
        if centroid_bounds.p_max[dim] == centroid_bounds.p_min[dim] {
            // Create leaf BVHBuildNode.
            new_leaf_node(arena, primitives, primitive_info, start, end, ordered_prims, bounds)
        } else {
            // Partition primitives based on split_method.
            let mid = match split_method {
                SplitMethod::Middle => {
                    let mut m = split_middle(primitive_info, start, end, dim, &centroid_bounds);

                    // For lots of prims with large overlapping bounding boxes, this may fail to partition; in that
                    // case use EqualCounts.
                    if m != start && m != end {
                        m = split_equal_counts(primitive_info, start, end, dim);
                    }
                    Some(m)
                }

                SplitMethod::EqualCounts => Some(split_equal_counts(primitive_info, start, end, dim)),

                SplitMethod::SAH => {
                    if n_primitives <= 2 {
                        // Partition primitives into equally-sized subsets.
                        Some(split_equal_counts(primitive_info, start, end, dim))
                    } else {
                        split_sah(
                            primitive_info,
                            start,
                            end,
                            n_primitives,
                            dim,
                            &centroid_bounds,
                            &bounds,
                            max_prims_in_node,
                        )
                    }
                }

                _ => panic!("sah::build(): Invalid split_method={:?}", split_method),
            };

            match mid {
                Some(mid) => {
                    // Create interior BVHBuildNode.
                    new_interior_node(
                        arena,
                        primitives,
                        split_method,
                        max_prims_in_node,
                        primitive_info,
                        start,
                        mid,
                        end,
                        dim,
                        total_nodes,
                        ordered_prims,
                    )
                }
                None => {
                    // Create leaf BVHBuildNode.
                    new_leaf_node(arena, primitives, primitive_info, start, end, ordered_prims, bounds)
                }
            }
        }
    }
}

/// Create interior `BVHBuildNode`.
///
/// * `arena`             - Shared arena for memory allocations.
/// * `primitives`        - The primitives in the node.
/// * `split_method`      - Middle|EqualCounts|SAH
/// * `max_prims_in_node` - Maximum number of primitives in the node.
/// * `primitive_info`    - Primitive information.
/// * `start`             - Starting index.
/// * `mid`               - Middle index.
/// * `end`               - Ending index.
/// * `dim`               - Split axis.
/// * `total_nodes`       - Used to return total number of nodes.
/// * `ordered_prims`     - Used to return a list of primitives ordered such that primitives in leaf nodes occupy
///                         contiguous ranges in the vector.
fn new_interior_node(
    arena: &SharedArena<BVHBuildNode>,
    primitives: &[ArcPrimitive],
    split_method: SplitMethod,
    max_prims_in_node: u8,
    primitive_info: &mut [BVHPrimitiveInfo],
    start: usize,
    mid: usize,
    end: usize,
    dim: Axis,
    total_nodes: &mut usize,
    ordered_prims: Arc<Mutex<Vec<ArcPrimitive>>>,
) -> ArenaArc<BVHBuildNode> {
    let c0 = build(
        arena,
        primitives,
        split_method,
        max_prims_in_node,
        primitive_info,
        start,
        mid,
        total_nodes,
        Arc::clone(&ordered_prims),
    );
    let c1 = build(
        arena,
        primitives,
        split_method,
        max_prims_in_node,
        primitive_info,
        mid,
        end,
        total_nodes,
        Arc::clone(&ordered_prims),
    );

    arena.alloc_arc(BVHBuildNode::new_interior_node(dim, c0, c1))
}

/// Create leaf `BVHBuildNode`.
///
/// * `arena`          - Shared arena for memory allocations.
/// * `primitives`     - The primitives in the node.
/// * `primitive_info` - Primitive information.
/// * `start`          - Starting index.
/// * `end`            - Ending index.
/// * `total_nodes`    - Used to return total number of nodes.
/// * `ordered_prims`  - Used to return a list of primitives ordered such that primitives in leaf nodes occupy
///                      contiguous ranges in the vector.
/// * `bounds`         - Bounds of all primitives in BVH node.
fn new_leaf_node(
    arena: &SharedArena<BVHBuildNode>,
    primitives: &[ArcPrimitive],
    primitive_info: &[BVHPrimitiveInfo],
    start: usize,
    end: usize,
    ordered_prims: Arc<Mutex<Vec<ArcPrimitive>>>,
    bounds: Bounds3f,
) -> ArenaArc<BVHBuildNode> {
    // Create leaf BVHBuildNode.
    let n_primitives = end - start;
    let prims = Arc::clone(&ordered_prims);
    let mut prims2 = prims.lock().expect("unable to lock ordered_prims");
    let first_prim_offset = prims2.len();
    for i in start..end {
        let prim_num = primitive_info[i].primitive_number;
        prims2.push(Arc::clone(&primitives[prim_num]));
    }

    arena.alloc_arc(BVHBuildNode::new_leaf_node(first_prim_offset, n_primitives, bounds))
}

/// Linear Bounding Volume Hierarchy using splitting planes that are midpoint of each region of space.
///
/// * `primitive_info`  - Vector containing all primitive info.
/// * `start`           - Starting index in primitive_info.
/// * `end`             - Ending index + 1 in primitive_info.
/// * `dim`             - Axis used to partition primitives.
/// * `centroid_bounds` - Bounding box of primtive centroids in primtive_info from start to end.
fn split_middle(
    primitive_info: &mut [BVHPrimitiveInfo],
    start: usize,
    end: usize,
    dim: Axis,
    centroid_bounds: &Bounds3f,
) -> usize {
    let pmid = (centroid_bounds.p_min[dim] + centroid_bounds.p_max[dim]) / 2.0;
    let infos = primitive_info[start..end].iter_mut();
    let split = itertools::partition(infos, |pi| pi.centroid[dim] < pmid);
    start + split
}

/// Partition primitives into equally sized subsets such that the first half of the primitives have smallest centroid
/// coordinate values along the chosen axis, and second have have the largest centroid coordinate values.
///
/// * `primitive_info`  - Vector containing all primitive info.
/// * `start`           - Starting index in primitive_info.
/// * `end`             - Ending index + 1 in primitive_info.
/// * `dim`             - Axis used to partition primitives.
fn split_equal_counts(primitive_info: &mut [BVHPrimitiveInfo], start: usize, end: usize, dim: Axis) -> usize {
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

/// Partition a subset of items between start and end inclusive such that:
/// - k^th element will be in its sorted order
/// - elements e in v[start, k - 1] will satisfy f(e, ek) == Ordering::Less
/// - elements e in v[k + 1, end] will satisfy f(e, ek) == Ordering::Greater and the k^th element is returned.
///
/// * `v`     - Vector to partition.
/// * `start` - Starting index.
/// * `end`   - Ending index + 1.
/// * `f`     - Predicate used for partitioning.
fn kth_element_by<F, T>(v: &mut [T], start: usize, k: usize, end: usize, f: F) -> Option<T>
where
    F: Fn(&T, &T) -> Ordering,
    T: Copy,
{
    if start >= end || end > v.len() || k < start || k >= end {
        return None;
    }

    let w = v[start..end].as_mut();
    let i = *kth_by(w, k - start, |&x, &y| f(&x, &y));

    Some(i)
}

/// Partition primitives using Surface Area Heuristic.
///
/// If the algorithm is able to partition primitives it will return the pivot index (mid) for interior node creation;
/// otherwise None is returned to indicate leaf node creation.
///
/// * `primitive_info`    - Vector containing all primitive info.
/// * `start`             - Start index in primitive_info.
/// * `end`               - End index in primitive_info.
/// * `n_primitives`      - Number of primitives between start and end.
/// * `dim`               - Axis used to partition primitives.
/// * `centroid_bounds`   - Bounding box of primtive centroids in primtive_info from start to end.
/// * `bounds`            - Bound box of all primitives in BVH node.
/// * `max_prims_in_node` - Maximum primitives allowed in node.
fn split_sah(
    primitive_info: &mut [BVHPrimitiveInfo],
    start: usize,
    end: usize,
    n_primitives: usize,
    dim: Axis,
    centroid_bounds: &Bounds3f,
    bounds: &Bounds3f,
    max_prims_in_node: u8,
) -> Option<usize> {
    // Partition primitives using approximate SAH.
    // Allocate BucketInfo for SAH partition buckets.
    let mut buckets = [BucketInfo::default(); N_BUCKETS];

    // Initialize BucketInfo for SAH partition buckets.
    for i in start..end {
        let mut b = (N_BUCKETS as Float * centroid_bounds.offset(&primitive_info[i].centroid)[dim]) as usize;

        if b == N_BUCKETS {
            b = N_BUCKETS - 1;
        }

        buckets[b].count += 1;
        buckets[b].bounds = buckets[b].bounds.union(&primitive_info[i].bounds);
    }

    // Compute costs for splitting after each bucket.
    let mut cost = [0.0_f32; N_BUCKETS - 1];
    for i in 0..N_BUCKETS - 1 {
        let (mut b0, mut b1) = (Bounds3f::EMPTY, Bounds3f::EMPTY);
        let (mut count0, mut count1) = (0, 0);

        for j in 0..=i {
            b0 = b0.union(&buckets[j].bounds);
            count0 += buckets[j].count;
        }

        for j in i + 1..N_BUCKETS {
            b1 = b1.union(&buckets[j].bounds);
            count1 += buckets[j].count;
        }

        cost[i] =
            1.0 + (count0 as Float * b0.surface_area() + count1 as Float * b1.surface_area()) / bounds.surface_area();
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

    // Either create leaf or split primitives at selected SAH bucket.
    let leaf_cost = n_primitives as Float;
    if n_primitives > max_prims_in_node as usize || min_cost < leaf_cost {
        // Partition primitives at selected SAH bucket and return the pivot point as mid.
        let infos = primitive_info[start..end].iter_mut();
        let split = itertools::partition(infos, |pi| {
            let mut b = (N_BUCKETS as Float * centroid_bounds.offset(&pi.centroid)[dim]) as usize;
            if b == N_BUCKETS {
                b = N_BUCKETS - 1;
            }
            b <= min_cost_split_bucket
        });
        let pmid = start + split;
        Some(pmid)
    } else {
        // No split occurred. Indicate creation of leaf node.
        None
    }
}
