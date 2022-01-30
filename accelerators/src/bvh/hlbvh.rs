//! HLBVH Algorithm

use super::common::*;
use super::morton::*;
use core::geometry::*;
use core::pbrt::*;
use core::primitive::*;
use rayon::prelude::*;
use shared_arena::{ArenaArc, SharedArena};
use std::cell::RefCell;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

const N_BUCKETS: usize = 12;
const FIRST_BIT_INDEX: usize = N_BITS - 1 - N_BUCKETS; // The index of the next bit to try splitting
const MORTON_BITS: u32 = 10;
const MORTON_SCALE: u32 = 1 << MORTON_BITS;

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
    let bounds = primitive_info
        .iter()
        .fold(Bounds3f::EMPTY, |b, pi| b.union(&pi.bounds));

    // Compute Morton indices of primitives.
    let morton_prims: Vec<MortonPrimitive> = primitive_info
        .par_iter()
        .map(|&pi| {
            let v = bounds.offset(&pi.centroid) * MORTON_SCALE as Float;
            let morton_code = encode_morton_3(&v);
            MortonPrimitive::new(pi.primitive_number, morton_code)
        })
        .collect();

    // Radix sort primitive Morton indices.
    let mut morton_prims_cell = RefCell::new(morton_prims);
    radix_sort(&mut morton_prims_cell);
    let morton_prims = morton_prims_cell.into_inner();

    // Create LBVH treelets at bottom of BVH.
    const MASK: u32 = 0b00111111111111000000000000000000;
    let mut treelets_to_build: Vec<(usize, usize)> = vec![]; // (start, n_primitives)
    let (mut start, mut end) = (0, 1);
    while end <= morton_prims.len() {
        if end == morton_prims.len()
            || ((morton_prims[start].morton_code & MASK) != (morton_prims[end].morton_code & MASK))
        {
            treelets_to_build.push((start, end - start));
            start = end;
        }

        end += 1;
    }

    // Create LBVHs for treelets in parallel.
    let atomic_total = AtomicUsize::new(0);
    let ordered_prims_offset = AtomicUsize::new(0);
    let mut treelets: Vec<ArenaArc<BVHBuildNode>> = treelets_to_build
        .par_iter()
        .map(|&(start_index, n_primitives)| {
            // Generate i^th LBVH treelet.
            let mut nodes_created = 0;
            let build_node = emit_lbvh(
                arena,
                primitives,
                max_prims_in_node as usize,
                primitive_info,
                &morton_prims[start_index..],
                n_primitives,
                &mut nodes_created,
                Arc::clone(&ordered_prims),
                &ordered_prims_offset,
                Some(FIRST_BIT_INDEX),
            );
            atomic_total.fetch_add(nodes_created, Ordering::SeqCst);
            build_node
        })
        .collect();

    *total_nodes = atomic_total.into_inner();

    // Create and return SAH BVH from LBVH treelets.
    build_upper_sah(
        arena,
        &mut treelets,
        0,
        treelets_to_build.len(),
        total_nodes,
    )
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
    bit_index: Option<usize>,
) -> ArenaArc<BVHBuildNode> {
    debug_assert!(n_primitives > 0);

    if bit_index.is_none() || n_primitives < max_prims_in_node {
        // Create and return leaf node of LBVH treelet.
        let mut bounds = Bounds3f::EMPTY;
        let first_prim_offset = ordered_prims_offset.fetch_add(n_primitives, Ordering::SeqCst);

        let prims = Arc::clone(&ordered_prims);
        let mut prims2 = prims.lock().expect("unabled to lock ordered_prims");
        for i in 0..n_primitives {
            let primitive_index = morton_prims[i].primitive_index;
            prims2[first_prim_offset + i] = Arc::clone(&primitives[primitive_index]);
            bounds = bounds.union(&primitive_info[primitive_index].bounds);
        }

        *total_nodes += 1;
        arena.alloc_arc(BVHBuildNode::new_leaf_node(
            first_prim_offset,
            n_primitives,
            bounds,
        ))
    } else if let Some(bit_idx) = bit_index {
        let mask = 1 << bit_idx;
        // Advance to next subtree level if there's no LBVH split for this bit
        if (morton_prims[0].morton_code & mask)
            == (morton_prims[n_primitives - 1].morton_code & mask)
        {
            return emit_lbvh(
                arena,
                primitives,
                max_prims_in_node,
                primitive_info,
                morton_prims,
                n_primitives,
                total_nodes,
                Arc::clone(&ordered_prims),
                ordered_prims_offset,
                Some(bit_idx - 1),
            );
        }

        // Find LBVH split point for this dimension
        let (mut search_start, mut search_end) = (0, n_primitives - 1);
        while search_start + 1 != search_end {
            debug_assert!(search_start != search_end);

            let mid = (search_start + search_end) / 2;
            if (morton_prims[search_start].morton_code & mask)
                == (morton_prims[mid].morton_code & mask)
            {
                search_start = mid;
            } else {
                debug_assert!(
                    morton_prims[mid].morton_code & mask
                        == morton_prims[search_end].morton_code & mask
                );
                search_end = mid;
            }
        }
        let split_offset = search_end;
        debug_assert!(split_offset < n_primitives - 1);
        debug_assert!(
            morton_prims[split_offset - 1].morton_code & mask
                != morton_prims[split_offset].morton_code & mask
        );

        // Create and return interior LBVH node
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
            Some(bit_idx - 1),
        );
        let c1 = emit_lbvh(
            arena,
            primitives,
            max_prims_in_node,
            primitive_info,
            &morton_prims[split_offset..],
            n_primitives - split_offset,
            total_nodes,
            Arc::clone(&ordered_prims),
            ordered_prims_offset,
            Some(bit_idx - 1),
        );

        arena.alloc_arc(BVHBuildNode::new_interior_node(
            Axis::from(bit_idx % 3),
            c0,
            c1,
        ))
    } else {
        panic!("HLBVH::emit_lbvh(): bit_index is none");
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
    debug_assert!(start < end);

    let n_nodes = end - start;
    if n_nodes == 1 {
        return ArenaArc::clone(&treelet_roots[start]);
    }
    *total_nodes += 1;

    // Compute bounds of all nodes under this HLBVH node
    let mut bounds = Bounds3f::EMPTY;
    for treelet_root in treelet_roots.iter().take(end).skip(start) {
        bounds = bounds.union(&treelet_root.bounds);
    }

    // Compute bound of HLBVH node centroids, choose split dimension dim.
    let mut centroid_bounds = Bounds3f::EMPTY;
    for treelet_root in treelet_roots.iter().take(end).skip(start) {
        let centroid = (treelet_root.bounds.p_min + treelet_root.bounds.p_max) * 0.5;
        centroid_bounds = centroid_bounds.union(&centroid);
    }

    let dim = centroid_bounds.maximum_extent();

    // Make sure the SAH split below does something... ?
    debug_assert!(centroid_bounds.p_max[dim] != centroid_bounds.p_min[dim]);

    // Allocate BucketInfo for SAH partition buckets
    let mut buckets = [BucketInfo::default(); N_BUCKETS];

    // Initialize BucketInfo for HLBVH SAH partition buckets
    for treelet_root in treelet_roots.iter().take(end).skip(start) {
        let centroid = (treelet_root.bounds.p_min[dim] + treelet_root.bounds.p_max[dim]) * 0.5;

        let mut b = (N_BUCKETS as Float
            * ((centroid - centroid_bounds.p_min[dim])
                / (centroid_bounds.p_max[dim] - centroid_bounds.p_min[dim])))
            as usize;
        if b == N_BUCKETS {
            b = N_BUCKETS - 1;
        }
        debug_assert!(b > 0);
        debug_assert!(b < N_BUCKETS);

        buckets[b].count += 1;
        buckets[b].bounds = buckets[b].bounds.union(&treelet_root.bounds);
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

        *cost_i = 0.125
            + (count0 as Float * b0.surface_area() + count1 as Float * b1.surface_area())
                / bounds.surface_area();
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
            * ((centroid - centroid_bounds.p_min[dim])
                / (centroid_bounds.p_max[dim] - centroid_bounds.p_min[dim])))
            as usize;
        if b == N_BUCKETS {
            b = N_BUCKETS - 1;
        }
        debug_assert!(b > 0);
        debug_assert!(b < N_BUCKETS);
        b <= min_cost_split_bucket
    });

    let mid = start + split;
    debug_assert!(mid > start);
    debug_assert!(mid < end);

    arena.alloc_arc(BVHBuildNode::new_interior_node(
        dim,
        build_upper_sah(arena, treelet_roots, start, mid, total_nodes),
        build_upper_sah(arena, treelet_roots, mid, end, total_nodes),
    ))
}
