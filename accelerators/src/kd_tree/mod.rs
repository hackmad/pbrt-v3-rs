//! KD Tree.

use core::geometry::*;
use core::interaction::*;
use core::light::*;
use core::material::*;
use core::paramset::*;
use core::pbrt::*;
use core::primitive::*;
use core::reflection::*;

mod common;
use common::*;

/// KD Tree Accelerator.
#[derive(Clone)]
pub struct KDTreeAccel {
    /// Intersection cost.
    pub isect_cost: i32,

    /// Traversal cost.
    pub traversal_cost: i32,

    /// Maximum number of primitives in leaf node.
    pub max_prims: usize,

    /// Bonus value used when one of the 2 regions along split are empty.
    pub empty_bonus: Float,

    /// The primitives.
    pub primitives: Vec<ArcPrimitive>,

    /// Indexes into `primitives`.
    pub primitive_indices: Vec<usize>,

    /// The tree root.
    pub nodes: Vec<KdAccelNode>,

    /// Next free node.
    pub next_free_node: usize,

    /// Bounding box.
    pub bounds: Bounds3f,
}

impl KDTreeAccel {
    /// Create a new KD Tree Accelerator.
    ///
    /// * `primitives`     - Pimitives.
    /// * `isect_cost`     - Intersection cost (default 80).
    /// * `traversal_cost` - Traversal cost (default 1).
    /// * `empty_bonus`    - Bonus value used when one of the 2 regions along
    ///                      split are empty (default 0.5).
    /// * `max_prims`      - Maximum number of primitives in leaf (default 1).
    /// * `max_depth`      - Maximum depth of tree (default -1).
    pub fn new(
        primitives: &[ArcPrimitive],
        isect_cost: i32,
        traversal_cost: i32,
        empty_bonus: Float,
        max_prims: usize,
        max_depth: i32,
    ) -> Self {
        // Build kd-tree for accelerator.
        let count = primitives.len();
        let next_free_node = 0;

        let max_depth = if max_depth <= 0 {
            (8.0 + 1.3 * (count as i64).log2int() as Float).round() as usize
        } else {
            max_depth as usize
        };

        // Compute bounds for kd-tree construction.
        let mut bounds = Bounds3f::EMPTY;
        let mut prim_bounds = Vec::<Bounds3f>::with_capacity(count);
        for prim in primitives.iter() {
            let b = prim.world_bound();
            bounds = bounds.union(&b);
            prim_bounds.push(b);
        }

        // Allocate working memory for kd-tree construction.
        let mut edges = [
            vec![BoundEdge::default(); 2 * count],
            vec![BoundEdge::default(); 2 * count],
            vec![BoundEdge::default(); 2 * count],
        ];

        let mut prims0: Vec<usize> = vec![0_usize; count];
        let mut prims1: Vec<usize> = vec![0_usize; (max_depth + 1) * count];

        // Initialize prim_nums for kd-tree construction.
        let prim_nums: Vec<usize> = (0..count).collect();

        // Start recursive construction of kd-tree
        let mut kd_tree = Self {
            isect_cost,
            traversal_cost,
            max_prims,
            empty_bonus,
            primitives: primitives.to_vec(),
            primitive_indices: vec![],
            nodes: vec![],
            next_free_node,
            bounds,
        };

        kd_tree.build_tree(
            0,
            &bounds,
            &prim_bounds,
            &prim_nums,
            max_depth,
            &mut edges,
            &mut prims0,
            &mut prims1,
            0,
        );

        kd_tree
    }

    /// Build the subtree.
    ///
    /// * `node_num`        - Node index.
    /// * `node_bounds`     - Bounds for the node.
    /// * `all_prim_bounds` - Bounds of all primitives.
    /// * `prim_nums`       - Primitive indices.
    /// * `depth`           - Depth (decreases from max_depth to 0 from root to leaf).
    /// * `edges`           - Projections of bounding boxes along split axes.
    ///                       NOTE: This could be statically allocated within the
    ///                       function but it is being passed in to reduce those
    ///                       allocations. Its use in each call does not affect
    ///                       other branches of recursion.
    /// * `prims0`          - Primitives below split axis.
    /// * `prims1`          - Primitives above split axis.
    /// * `bad_refines`     - Keeps track of how many bad splits have been made
    ///                       so far above the current node of the tree.
    fn build_tree(
        &mut self,
        node_num: usize,
        node_bounds: &Bounds3f,
        all_prim_bounds: &[Bounds3f],
        prim_nums: &[usize],
        depth: usize,
        edges: &mut [Vec<BoundEdge>; 3],
        prims0: &mut [usize],
        prims1: &mut [usize],
        bad_refines: u8,
    ) {
        assert_eq!(node_num, self.next_free_node);

        let mut bad_refines: u8 = bad_refines;

        let n_primitives = prim_nums.len();

        // Get next free node from nodes.
        let n_alloced_noded = self.nodes.len();
        if self.next_free_node == n_alloced_noded {
            let n = max(2 * n_alloced_noded, 512);
            self.nodes.resize_with(n, KdAccelNode::default);
        }
        self.next_free_node += 1;

        // Initialize leaf node if termination criteria met.
        if n_primitives <= self.max_prims || depth == 0 {
            self.nodes[node_num] = KdAccelNode::new_leaf(prim_nums, &mut self.primitive_indices);
            return;
        }

        // Initialize interior node and continue recursion.

        // Choose split axis position for interior node.
        let mut best_axis: Option<u8> = None;
        let mut best_offset: isize = -1;
        let mut best_cost = INFINITY;
        let old_cost = self.isect_cost as Float * n_primitives as Float;
        let total_sa = node_bounds.surface_area();
        let inv_total_sa = 1.0 / total_sa;
        let d = node_bounds.diagonal();

        // Choose which axis to split.
        let mut axis = node_bounds.maximum_extent() as usize;
        let mut retries = 0_u8;
        loop {
            // Initialize edges for `axis`.
            for (i, pn) in prim_nums.iter().enumerate().take(n_primitives) {
                let bounds = &all_prim_bounds[*pn];
                edges[axis][2 * i] = BoundEdge::new(bounds.p_min[axis], *pn, EdgeType::Start);
                edges[axis][2 * i + 1] = BoundEdge::new(bounds.p_max[axis], *pn, EdgeType::End);
            }

            // Sort `edges` for `axis`.
            edges[axis][0..2 * n_primitives].sort_by(|e0, e1| e0.partial_cmp(e1).unwrap());

            // Compute cost of all splits for `axis` to find best.
            let (mut n_below, mut n_above) = (0, n_primitives);
            for i in 0..2 * n_primitives {
                if edges[axis][i].edge_type == EdgeType::End {
                    n_above -= 1;
                }

                let edge_t = edges[axis][i].t;
                if edge_t > node_bounds.p_min[axis] && edge_t < node_bounds.p_max[axis] {
                    // Compute cost for split at i^th edge.

                    // Compute child surface areas for split at `edge_t`.
                    let other_axis_0 = (axis + 1) % 3;
                    let other_axis_1 = (axis + 2) % 3;
                    let below_sa = 2.0
                        * (d[other_axis_0] * d[other_axis_1]
                            + (edge_t - node_bounds.p_min[axis])
                                * (d[other_axis_0] + d[other_axis_1]));
                    let above_sa = 2.0
                        * (d[other_axis_0] * d[other_axis_1]
                            + (node_bounds.p_max[axis] - edge_t)
                                * (d[other_axis_0] + d[other_axis_1]));
                    let p_below = below_sa * inv_total_sa;
                    let p_above = above_sa * inv_total_sa;
                    let eb = if n_above == 0 || n_below == 0 {
                        self.empty_bonus
                    } else {
                        0.0
                    };
                    let cost = self.traversal_cost as Float
                        + (self.isect_cost as Float)
                            * (1.0 - eb)
                            * (p_below * n_below as Float + p_above * n_above as Float);

                    // Update best split if this is lowest cost so far
                    if cost < best_cost {
                        best_cost = cost;
                        best_axis = Some(axis as u8);
                        best_offset = i as isize;
                    }
                }
                if edges[axis][i].edge_type == EdgeType::Start {
                    n_below += 1;
                }
            }
            assert!(n_below == n_primitives && n_above == 0);

            // Create leaf node if no good splits were found.
            if best_axis.is_none() && retries < 2 {
                retries += 1;
                axis = (axis + 1) % 3;
                // retry split
            } else {
                break;
            }
        }

        if best_cost > old_cost {
            bad_refines += 1
        }

        if (best_cost > 4.0 * old_cost && n_primitives < 16)
            || best_axis.is_none()
            || bad_refines == 3
        {
            self.nodes[node_num] = KdAccelNode::new_leaf(prim_nums, &mut self.primitive_indices);
            return;
        }

        // Classify primitives with respect to split.
        let (mut n0, mut n1) = (0_usize, 0_usize);
        let best_axis = best_axis.unwrap() as usize;
        let best_offset = best_offset as usize;

        for i in 0..best_offset {
            if edges[best_axis][i].edge_type == EdgeType::Start {
                prims0[n0] = edges[best_axis][i].prim_num;
                n0 += 1;
            }
        }

        for i in (best_offset + 1)..(2 * n_primitives) {
            if edges[best_axis][i].edge_type == EdgeType::End {
                prims1[n1] = edges[best_axis][i].prim_num;
                n1 += 1;
            }
        }

        // Recursively initialize children nodes
        let t_split = edges[best_axis][best_offset].t;

        let mut bounds0 = *node_bounds;
        let mut bounds1 = *node_bounds;
        bounds0.p_max[best_axis] = t_split;
        bounds1.p_min[best_axis] = t_split;

        let prim_nums = prims0[0..n0].to_vec();
        self.build_tree(
            node_num + 1,
            &bounds0,
            all_prim_bounds,
            &prim_nums,
            depth - 1,
            edges,
            prims0,
            &mut prims1[n_primitives..],
            bad_refines,
        );

        let above_child = self.next_free_node;
        self.nodes[node_num] = KdAccelNode::new_interior(best_axis.into(), above_child, t_split);

        let prim_nums = prims1[0..n1].to_vec();
        self.build_tree(
            above_child,
            &bounds1,
            all_prim_bounds,
            &prim_nums,
            depth - 1,
            edges,
            prims0,
            &mut prims1[n_primitives..],
            bad_refines,
        );
    }
}

/// Tag `KDTreeAccel` as an `Aggregate`.
impl Aggregate for KDTreeAccel {}

impl Primitive for KDTreeAccel {
    /// Returns a bounding box in the world space.
    fn world_bound(&self) -> Bounds3f {
        self.bounds
    }

    /// Returns geometric details if a ray intersects the primitive and updates
    /// the t_max parameter of the ray. If there is no intersection, `None` is
    /// returned.
    ///
    /// * `r`                  - The ray.
    fn intersect(&self, r: &mut Ray) -> Option<SurfaceInteraction> {
        if self.nodes.is_empty() {
            return None;
        }

        // Compute initial parametric range of ray inside kd-tree extent.
        let bounds_hit = self.bounds.intersect_p(r);
        if bounds_hit.is_none() {
            return None;
        }
        let (mut t_min, mut t_max) = bounds_hit.unwrap();

        // Prepare to traverse kd-tree for ray.
        let inv_dir = Vector3f::new(1.0 / r.d.x, 1.0 / r.d.y, 1.0 / r.d.z);

        let mut todo = [KdToDo::default(); MAX_TO_DO];
        let mut todo_pos = 0;
        let mut node_idx = 0;
        let mut node = self.nodes.get(node_idx);

        // Traverse kd-tree nodes in order for ray.
        let mut si: Option<SurfaceInteraction> = None;
        while node.is_some() {
            // Bail out if we found a hit closer than the current node.
            if r.t_max < t_min {
                break;
            }

            match node {
                Some(&KdAccelNode::Interior {
                    split_pos,
                    split_axis,
                    above_child,
                }) => {
                    // Process kd-tree interior node.
                    // Compute parametric distance along ray to split plane.
                    let t_plane = (split_pos - r.o[split_axis]) * inv_dir[split_axis];

                    // Get node children pointers for ray.
                    let below_first = (r.o[split_axis] < split_pos)
                        || (r.o[split_axis] == split_pos && r.d[split_axis] <= 0.0);

                    let (first_child_idx, second_child_idx) = if below_first {
                        (node_idx + 1, above_child)
                    } else {
                        (above_child, node_idx + 1)
                    };
                    let first_child = self.nodes.get(first_child_idx);
                    let second_child = self.nodes.get(second_child_idx);

                    // Advance to next child node, possibly enqueue other child.
                    if t_plane > t_max || t_plane <= 0.0 {
                        node = first_child;
                        node_idx = first_child_idx;
                    } else if t_plane < t_min {
                        node = second_child;
                        node_idx = second_child_idx;
                    } else {
                        // Enqueue `second_child` in todo list.
                        todo[todo_pos].node = second_child;
                        todo[todo_pos].node_idx = second_child_idx;
                        todo[todo_pos].t_min = t_plane;
                        todo[todo_pos].t_max = t_max;
                        todo_pos += 1;

                        node = first_child;
                        node_idx = first_child_idx;
                        t_max = t_plane;
                    }
                }

                Some(&KdAccelNode::Leaf {
                    n_prims,
                    one_primitive,
                    primitive_indices_offset,
                }) => {
                    // Check for intersections inside leaf node.
                    if n_prims == 1 {
                        // Check one primitive inside leaf node.
                        if let Some(si2) = self.primitives[one_primitive].intersect(r) {
                            si = Some(si2);
                        }
                    } else {
                        for i in 0..n_prims {
                            // Check one primitive inside leaf node.
                            let index = self.primitive_indices[primitive_indices_offset + i];
                            if let Some(si2) = self.primitives[index].intersect(r) {
                                si = Some(si2);
                            }
                        }
                    }

                    // Grab next node to process from todo list.
                    if todo_pos > 0 {
                        todo_pos -= 1;
                        KdToDo {
                            node,
                            node_idx,
                            t_min,
                            t_max,
                        } = todo[todo_pos];
                    } else {
                        break;
                    }
                }

                None => {
                    break;
                }
            }
        } // while

        si
    }

    /// Returns `true` if a ray-primitive intersection succeeds; otherwise `false`.
    ///
    /// * `r`                  - The ray.
    fn intersect_p(&self, r: &Ray) -> bool {
        if self.nodes.is_empty() {
            return false;
        }

        // Compute initial parametric range of ray inside kd-tree extent.
        let bounds_hit = self.bounds.intersect_p(r);
        if bounds_hit.is_none() {
            return false;
        }
        let (mut t_min, mut t_max) = bounds_hit.unwrap();

        // Prepare to traverse kd-tree for ray.
        let inv_dir = Vector3f::new(1.0 / r.d.x, 1.0 / r.d.y, 1.0 / r.d.z);
        let mut todo = [KdToDo::default(); MAX_TO_DO];
        let mut todo_pos = 0;
        let mut node_idx = 0;
        let mut node = self.nodes.get(node_idx);

        // Traverse kd-tree nodes in order for ray.
        while node.is_some() {
            // Bail out if we found a hit closer than the current node.
            if r.t_max < t_min {
                break;
            }

            match node {
                Some(&KdAccelNode::Interior {
                    split_pos,
                    split_axis,
                    above_child,
                }) => {
                    // Process kd-tree interior node.
                    // Compute parametric distance along ray to split plane.
                    let axis = split_axis as usize;
                    let t_plane = (split_pos - r.o[axis]) * inv_dir[axis];

                    // Get node children pointers for ray.
                    let below_first =
                        (r.o[axis] < split_pos) || (r.o[axis] == split_pos && r.d[axis] <= 0.0);

                    let (first_child_idx, second_child_idx) = if below_first {
                        (node_idx + 1, above_child)
                    } else {
                        (above_child, node_idx + 1)
                    };
                    let first_child = self.nodes.get(first_child_idx);
                    let second_child = self.nodes.get(second_child_idx);

                    // Advance to next child node, possibly enqueue other child.
                    if t_plane > t_max || t_plane <= 0.0 {
                        node = first_child;
                        node_idx = first_child_idx;
                    } else if t_plane < t_min {
                        node = second_child;
                        node_idx = second_child_idx;
                    } else {
                        // Enqueue `second_child` in todo list.
                        todo[todo_pos].node = second_child;
                        todo[todo_pos].node_idx = second_child_idx;
                        todo[todo_pos].t_min = t_plane;
                        todo[todo_pos].t_max = t_max;
                        todo_pos += 1;

                        node = first_child;
                        node_idx = first_child_idx;
                        t_max = t_plane;
                    }
                }
                Some(&KdAccelNode::Leaf {
                    n_prims,
                    one_primitive,
                    primitive_indices_offset,
                }) => {
                    // Check for intersections inside leaf node.
                    if n_prims == 1 {
                        // Check one primitive inside leaf node.
                        if self.primitives[one_primitive].intersect_p(r) {
                            return true;
                        }
                    } else {
                        for i in 0..n_prims {
                            // Check one primitive inside leaf node.
                            let index = self.primitive_indices[primitive_indices_offset + i];
                            if self.primitives[index].intersect_p(r) {
                                return true;
                            }
                        }
                    }

                    // Grab next node to process from todo list.
                    if todo_pos > 0 {
                        todo_pos -= 1;
                        KdToDo {
                            node,
                            node_idx,
                            t_min,
                            t_max,
                        } = todo[todo_pos];
                    } else {
                        break;
                    }
                }

                None => {
                    break;
                }
            }
        }

        false
    }

    /// Returns a reference to the AreaLight that describes the primitive’s
    /// emission distribution, if the primitive is itself a light source.
    /// If the primitive is not emissive, this method should return `None`.
    ///
    /// *NOTE*: This should never be called. Calling code should directly call
    /// get_area_light() on the primitive from the ray-primitive intersection.
    fn get_area_light(&self) -> Option<ArcLight> {
        error!(
            "TransformedPrimitive::get_area_light() shouldn't be called; \
            should've gone to GeometricPrimitive."
        );
        None
    }

    /// Returns a reference to the material instance assigned to the primitive.
    /// If `None` is returned, ray intersections with the primitive should be
    /// ignored; the primitive only serves to delineate a volume of space for
    /// participating media. This method is also used to check if two rays have
    /// intersected the same object by comparing their Material pointers.
    ///
    /// *NOTE*: This should never be called. Calling code should directly call
    /// get_material() on the primitive from the ray-primitive intersection.
    fn get_material(&self) -> Option<ArcMaterial> {
        error!(
            "TransformedPrimitive::get_material() shouldn't be called; \
            should've gone to GeometricPrimitive."
        );
        None
    }

    /// Initializes representations of the light-scattering properties of the
    /// material at the intersection point on the surface.
    ///
    /// *NOTE*: This should never be called. Calling code should directly call
    /// compute_scattering_functions() on the primitive from the ray-primitive
    /// intersection.
    ///
    /// * `_si`                   - The surface interaction at the intersection.
    /// * `_mode`                 - Transport mode.
    /// * `_allow_multiple_lobes` - Allow multiple lobes.
    /// * `_bsdf`                 - The computed BSDF.
    /// * `_bssrdf`               - The computed BSSSRDF.
    fn compute_scattering_functions<'scene>(
        &self,
        _si: &mut SurfaceInteraction<'scene>,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
        _bsdf: &mut Option<BSDF>,
        _bssrdf: &mut Option<BSDF>,
    ) {
        error!(
            "TransformedPrimitive::compute_scattering_functions() shouldn't be \
            called; should've gone to GeometricPrimitive."
        );
    }
}

impl From<(&ParamSet, &[ArcPrimitive])> for KDTreeAccel {
    /// Create a `KDTreeAccel ` from given parameter set and primitives.
    ///
    /// * `p` - Tuple containing the parameter set and primitives.
    fn from(p: (&ParamSet, &[ArcPrimitive])) -> Self {
        let (params, prims) = p;
        let isect_cost = params.find_one_int("intersectcost", 80);
        let trav_cost = params.find_one_int("traversalcost", 1);
        let empty_bonus = params.find_one_float("emptybonus", 0.5);
        let max_prims = params.find_one_int("maxprims", 1) as usize;
        let max_depth = params.find_one_int("maxdepth", -1);

        Self::new(
            prims,
            isect_cost,
            trav_cost,
            empty_bonus,
            max_prims,
            max_depth,
        )
    }
}
