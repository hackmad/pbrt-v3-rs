//! HLBVH Algorithm

use super::common::*;
use super::morton::*;
use core::app::OPTIONS;
use core::geometry::*;
use core::pbrt::*;
use core::primitive::*;
use itertools::Itertools;
use shared_arena::{ArenaArc, SharedArena};
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

const N_BUCKETS: usize = 12;
const FIRST_BIT_INDEX: usize = N_BITS - 1 - N_BUCKETS; // The index of the next bit to try splitting
const MORTON_BITS: u32 = 10;
const MORTON_SCALE: Float = (1 << MORTON_BITS) as Float;

/// Build the BVH structure using HLBVH algorithm.
///
/// * `arena`             - Shared arena for memory allocations.
/// * `primitives`        - The primitives in the node.
/// * `max_prims_in_node` - Maximum number of primitives in the node.
/// * `primitive_info`    - Primitive information.
/// * `total_nodes`       - Used to return total number of nodes.
/// * `ordered_prims`     - Used to return a list of primitives ordered such that
///                         primitives in leaf nodes occupy contiguous ranges in
///                         the vector.
pub fn build(
    arena: &SharedArena<BVHBuildNode>,
    primitives: &[ArcPrimitive],
    max_prims_in_node: u8,
    primitive_info: &mut Vec<BVHPrimitiveInfo>,
    total_nodes: &mut usize,
    ordered_prims: Arc<Mutex<Vec<ArcPrimitive>>>,
) -> ArenaArc<BVHBuildNode> {
    // Compute bounds of all primitives in BVH node.
    let bounds = primitive_info.iter().fold(Bounds3f::EMPTY, |b, pi| b.union(&pi.bounds));

    // Compute Morton indices of primitives.
    let morton_prims = compute_morton_primitives(primitive_info, &bounds);

    // Radix sort primitive Morton indices.
    let mut morton_prims_cell = RefCell::new(morton_prims);
    radix_sort(&mut morton_prims_cell);
    let morton_prims = morton_prims_cell.into_inner();

    // Create LBVH treelets at bottom of BVH.
    // Find intervals of primitives for each treelet.
    const MASK: u32 = 0b00111111111111000000000000000000;
    let mut treelets_to_build: Vec<LBVHTreelet> = vec![];
    let mut start = 0_usize;
    let mut end = 1_usize;
    while end <= morton_prims.len() {
        if end == morton_prims.len()
            || ((morton_prims[start].morton_code & MASK) != (morton_prims[end].morton_code & MASK))
        {
            // Add entry to `treelets_to_build` for this treelet.
            let n_primitives = end - start;
            treelets_to_build.push(LBVHTreelet::new(start, n_primitives));
            start = end;
        }

        end += 1;
    }

    // Create LBVHs for treelets in parallel.
    let atomic_total = AtomicUsize::new(0);
    let ordered_prims_offset = AtomicUsize::new(0);

    // Need to pre-allocate ordered_prims so just clone reference to first
    // primitive to use as default.
    {
        let mut prims = ordered_prims.lock().expect("Unable to lock ordered_prims");
        prims.resize_with(primitives.len(), || Arc::clone(&primitives[0]));
    }

    let mut treelets = build_treelets(
        arena,
        &treelets_to_build,
        primitives,
        max_prims_in_node,
        primitive_info,
        &morton_prims,
        ordered_prims,
        &ordered_prims_offset,
        &atomic_total,
    );

    *total_nodes = atomic_total.into_inner();

    // Create and return SAH BVH from LBVH treelets.
    build_upper_sah(arena, &mut treelets, 0, treelets_to_build.len(), total_nodes)
}

/// Compute Morton indices of primitives.
///
/// * `primitive_info` - Primitive information.
/// * `bounds`         - Bounds of all primitives in BVH node.
fn compute_morton_primitives(primitive_info: &Vec<BVHPrimitiveInfo>, bounds: &Bounds3f) -> Vec<MortonPrimitive> {
    let n = primitive_info.len();
    let morton_prims = Arc::new(Mutex::new(vec![MortonPrimitive::default(); n]));

    let n_threads = OPTIONS.threads();

    thread::scope(|scope| {
        let (tx, rx) = crossbeam_channel::bounded::<(usize, &BVHPrimitiveInfo)>(n_threads);

        for _ in 0..n_threads {
            let rxc = rx.clone();
            let morton_prims = Arc::clone(&morton_prims);
            scope.spawn(move || {
                for (i, pi) in rxc.iter() {
                    let centroid_offset = bounds.offset(&pi.centroid);
                    let v = centroid_offset * MORTON_SCALE;

                    let mut mp = morton_prims.lock().unwrap();
                    (*mp)[i].primitive_index = pi.primitive_number;
                    (*mp)[i].morton_code = encode_morton_3(&v);
                }
            });
        }
        drop(rx); // Drop extra rx since we've cloned one for each worker.

        // Send work.
        for (i, pi) in primitive_info.iter().enumerate() {
            tx.send((i, pi)).unwrap();
        }
    });

    let mut mp = morton_prims.lock().unwrap();
    std::mem::replace(&mut mp, vec![])
}

/// Create LBVHs for treelets in parallel.
///
/// * `arena`                - Shared arena for memory allocations.
/// * `treelets_to_build`    - Treelets to build.
/// * `primitives`           - The primitives in the node.
/// * `max_prims_in_node`    - Maximum number of primitives in the node.
/// * `primitive_info`       - Primitive information.
/// * `morton_prims`         - Morton primitives.
/// * `ordered_prims`        - Used to return a list of primitives ordered such that
///                            primitives in leaf nodes occupy contiguous ranges in
///                            the vector.
/// * `ordered_prims_offset  - Index in `ordered_prims` for start of this node.
/// * `atomic_total`         - Used to return total of nodes created.
fn build_treelets(
    arena: &SharedArena<BVHBuildNode>,
    treelets_to_build: &[LBVHTreelet],
    primitives: &[ArcPrimitive],
    max_prims_in_node: u8,
    primitive_info: &[BVHPrimitiveInfo],
    morton_prims: &[MortonPrimitive],
    ordered_prims: Arc<Mutex<Vec<ArcPrimitive>>>,
    ordered_prims_offset: &AtomicUsize,
    atomic_total: &AtomicUsize,
) -> Vec<ArenaArc<BVHBuildNode>> {
    let n = treelets_to_build.len();
    let treelets: Arc<Mutex<BTreeMap<usize, ArenaArc<BVHBuildNode>>>> = Arc::new(Mutex::new(BTreeMap::new()));

    let n_threads = OPTIONS.threads();

    thread::scope(|scope| {
        let (tx, rx) = crossbeam_channel::bounded(n_threads);

        for _ in 0..n_threads {
            let rxc = rx.clone();
            let treelets = Arc::clone(&treelets);
            let ordered_prims = &ordered_prims;
            scope.spawn(move || {
                for i in rxc.iter() {
                    let tr: &LBVHTreelet = &treelets_to_build[i];

                    // Generate i^th LBVH treelet.
                    let mut nodes_created = 0_usize;
                    let end_index = tr.start_index + tr.n_primitives;
                    let build_node = emit_lbvh(
                        arena,
                        primitives,
                        max_prims_in_node as usize,
                        primitive_info,
                        &morton_prims[tr.start_index..end_index],
                        tr.n_primitives,
                        &mut nodes_created,
                        Arc::clone(ordered_prims),
                        &ordered_prims_offset,
                        FIRST_BIT_INDEX as isize,
                    );
                    atomic_total.fetch_add(nodes_created, Ordering::AcqRel);

                    let mut tl = treelets.lock().unwrap();
                    (*tl).insert(i, build_node);
                }
            });
        }
        drop(rx); // Drop extra rx since we've cloned one for each worker.

        // Send work.
        for i in 0..n {
            tx.send(i).unwrap();
        }
    });

    let tl = treelets.lock().unwrap();
    tl.values().map(ArenaArc::clone).collect_vec()
}

/// Builds a treelet by taking primitives with centroids in some region of space
/// and successively partitions them with splitting planes that divide the current
/// region of space into two halves along the center of the region along one of
/// the three axes.
///
/// * `arena`                - Shared arena for memory allocations.
/// * `primitives`           - The primitives.
/// * `max_prims_in_node`    - Maximum number of primitives in the node.
/// * `primitive_info`       - Primitive information.
/// * `morton_prims`         - Morton codes for primitives.
/// * `n_primitives`         - Number of primitives.
/// * `total_nodes`          - Total number of nodes.
/// * `ordered_prims`        - Used to return a list of primitives ordered such that
///                            primitives in leaf nodes occupy contiguous ranges in
///                            the vector.
/// * `ordered_prims_offset  - Index in `ordered_prims` for start of this node.
/// * `bit_index`            - The bit index.
fn emit_lbvh(
    arena: &SharedArena<BVHBuildNode>,
    primitives: &[ArcPrimitive],
    max_prims_in_node: usize,
    primitive_info: &[BVHPrimitiveInfo],
    morton_prims: &[MortonPrimitive],
    n_primitives: usize,
    total_nodes: &mut usize,
    ordered_prims: Arc<Mutex<Vec<ArcPrimitive>>>,
    ordered_prims_offset: &AtomicUsize,
    bit_index: isize,
) -> ArenaArc<BVHBuildNode> {
    assert!(n_primitives > 0);

    if bit_index == -1 || n_primitives < max_prims_in_node {
        // Create and return leaf node of LBVH treelet.
        let mut bounds = Bounds3f::EMPTY;
        let first_prim_offset = ordered_prims_offset.fetch_add(n_primitives, Ordering::AcqRel);

        let mut prims = ordered_prims.lock().expect("unabled to lock ordered_prims");
        for i in 0..n_primitives {
            let primitive_index = morton_prims[i].primitive_index;
            prims[first_prim_offset + i] = Arc::clone(&primitives[primitive_index]);
            bounds = bounds.union(&primitive_info[primitive_index].bounds);
        }

        *total_nodes += 1;
        arena.alloc_arc(BVHBuildNode::new_leaf_node(first_prim_offset, n_primitives, bounds))
    } else {
        let mask = 1 << bit_index;
        // Advance to next subtree level if there's no LBVH split for this bit.
        if (morton_prims[0].morton_code & mask) == (morton_prims[n_primitives - 1].morton_code & mask) {
            return emit_lbvh(
                arena,
                primitives,
                max_prims_in_node,
                primitive_info,
                morton_prims,
                n_primitives,
                total_nodes,
                ordered_prims,
                ordered_prims_offset,
                bit_index - 1,
            );
        }

        // Find LBVH split point for this dimension
        let mut search_start = 0_usize;
        let mut search_end = n_primitives - 1;
        while search_start + 1 != search_end {
            assert_ne!(search_start, search_end);

            let mid = (search_start + search_end) / 2;
            if (morton_prims[search_start].morton_code & mask) == (morton_prims[mid].morton_code & mask) {
                search_start = mid;
            } else {
                assert_eq!(
                    morton_prims[mid].morton_code & mask,
                    morton_prims[search_end].morton_code & mask
                );
                search_end = mid;
            }
        }
        let split_offset = search_end;
        assert!(split_offset <= n_primitives - 1);
        assert_ne!(
            morton_prims[split_offset - 1].morton_code & mask,
            morton_prims[split_offset].morton_code & mask
        );

        // Create and return interior LBVH node.
        let c0 = emit_lbvh(
            arena,
            primitives,
            max_prims_in_node,
            primitive_info,
            morton_prims,
            split_offset,
            total_nodes,
            Arc::clone(&ordered_prims),
            ordered_prims_offset,
            bit_index - 1,
        );
        let c1 = emit_lbvh(
            arena,
            primitives,
            max_prims_in_node,
            primitive_info,
            &morton_prims[split_offset..],
            n_primitives - split_offset,
            total_nodes,
            ordered_prims,
            ordered_prims_offset,
            bit_index - 1,
        );

        *total_nodes += 1;
        arena.alloc_arc(BVHBuildNode::new_interior_node(
            Axis::from((bit_index as usize) % 3),
            c0,
            c1,
        ))
    }
}

/// Creates a BVH of all the treelets.
///
/// * `arena`         - Shared arena for memory allocations.
/// * `treelet_roots` - Treelet roots.
/// * `start`         - Starting index. For first call it should be 0.
/// * `end`           - Ending index + 1. For first call it should be # of nodes.
/// * `total_nodes`   - Total number of nodes.
fn build_upper_sah(
    arena: &SharedArena<BVHBuildNode>,
    treelet_roots: &mut Vec<ArenaArc<BVHBuildNode>>,
    start: usize,
    end: usize,
    total_nodes: &mut usize,
) -> ArenaArc<BVHBuildNode> {
    assert!(start < end);

    let n_nodes = end - start;
    if n_nodes == 1 {
        return ArenaArc::clone(&treelet_roots[start]);
    }

    // Compute bounds of all nodes under this HLBVH node
    let mut bounds = Bounds3f::EMPTY;
    for i in start..end {
        bounds = bounds.union(&treelet_roots[i].bounds);
    }

    // Compute bound of HLBVH node centroids, choose split dimension `dim`.
    let mut centroid_bounds = Bounds3f::EMPTY;
    for i in start..end {
        let centroid = (treelet_roots[i].bounds.p_min + treelet_roots[i].bounds.p_max) * 0.5;
        centroid_bounds = centroid_bounds.union(&centroid);
    }

    let dim = centroid_bounds.maximum_extent();

    // FIXME: if this hits, what do we need to do?
    // Make sure the SAH split below does something... ?
    assert_ne!(centroid_bounds.p_max[dim], centroid_bounds.p_min[dim]);

    // Allocate BucketInfo for SAH partition buckets
    let mut buckets = [BucketInfo::default(); N_BUCKETS];

    // Initialize BucketInfo for HLBVH SAH partition buckets
    for i in start..end {
        let centroid = (treelet_roots[i].bounds.p_min[dim] + treelet_roots[i].bounds.p_max[dim]) * 0.5;

        let mut b = (N_BUCKETS as Float
            * ((centroid - centroid_bounds.p_min[dim]) / (centroid_bounds.p_max[dim] - centroid_bounds.p_min[dim])))
            as usize;
        if b == N_BUCKETS {
            b = N_BUCKETS - 1;
        }
        assert!(b < N_BUCKETS);

        buckets[b].count += 1;
        buckets[b].bounds = buckets[b].bounds.union(&treelet_roots[i].bounds);
    }

    // Compute costs for splitting after each bucket
    let mut cost = [0.0; N_BUCKETS - 1];
    for (i, cost_i) in cost.iter_mut().enumerate().take(N_BUCKETS - 1) {
        let (mut b0, mut b1) = (Bounds3f::EMPTY, Bounds3f::EMPTY);
        let (mut count0, mut count1) = (0, 0);

        for bucket in buckets.iter().take(i + 1) {
            b0 = b0.union(&bucket.bounds);
            count0 += bucket.count;
        }

        for bucket in buckets.iter().take(N_BUCKETS).skip(i + 1) {
            b1 = b1.union(&bucket.bounds);
            count1 += bucket.count;
        }

        *cost_i =
            0.125 + (count0 as Float * b0.surface_area() + count1 as Float * b1.surface_area()) / bounds.surface_area();
    }

    // Find bucket to split at that minimizes SAH metric
    let mut min_cost = cost[0];
    let mut min_cost_split_bucket = 0;
    for (i, cost_i) in cost.iter().enumerate().take(N_BUCKETS - 1).skip(1) {
        if *cost_i < min_cost {
            min_cost = *cost_i;
            min_cost_split_bucket = i;
        }
    }

    // Split nodes and create interior HLBVH SAH node
    let roots = treelet_roots[start..end].iter_mut();
    let split = itertools::partition(roots, |node| {
        let centroid = (node.bounds.p_min[dim] + node.bounds.p_max[dim]) * 0.5;
        let mut b = (N_BUCKETS as Float
            * ((centroid - centroid_bounds.p_min[dim]) / (centroid_bounds.p_max[dim] - centroid_bounds.p_min[dim])))
            as usize;
        if b == N_BUCKETS {
            b = N_BUCKETS - 1;
        }
        assert!(b < N_BUCKETS);
        b <= min_cost_split_bucket
    });

    let mid = start + split;
    assert!(mid > start);
    assert!(mid < end);

    *total_nodes += 1;
    arena.alloc_arc(BVHBuildNode::new_interior_node(
        dim,
        build_upper_sah(arena, treelet_roots, start, mid, total_nodes),
        build_upper_sah(arena, treelet_roots, mid, end, total_nodes),
    ))
}

/// Treelet to build.
struct LBVHTreelet {
    /// Starting index.
    start_index: usize,

    /// Number of primitives.
    n_primitives: usize,
}

impl LBVHTreelet {
    /// Create a new instance of `LBVHTreelet`.
    ///
    /// * `start_index`  - Starting index.
    /// * `n_primitives` - Number of primitives.
    fn new(start_index: usize, n_primitives: usize) -> Self {
        Self {
            start_index,
            n_primitives,
        }
    }
}
