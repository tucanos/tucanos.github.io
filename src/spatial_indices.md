# Spatial indices

`tmesh` provides two complementary spatial indices for search, projection, and interpolation workflows.

## Point index

Point index uses a kd-tree over mesh vertices.

- library: [kdtree](https://crates.io/crates/kdtree),
- nearest vertex queries,
- ball/radius queries.

## Element/object index

Object index uses a BVH over element bounding boxes.

- library: [parry](https://parry.rs/) ([parry2d](https://crates.io/crates/parry2d-f64), [parry3d](https://crates.io/crates/parry3d-f64)),
- nearest element queries,
- closest-point projections,
- projected distance and projected coordinates.

This is the core index used for barycentric point location and linear interpolation.

## Interpolation methods

Two mesh-field interpolation strategies are used.

Nearest-neighbor interpolation:

- uses point index,
- assigns each query point the value at its nearest source vertex.

Linear simplex interpolation:

- uses object index to select a nearby source simplex,
- computes barycentric coordinates $\lambda_i$,
- interpolates with

$$
f(\mathbf{x}) = \sum_{i=1}^{n_v} \lambda_i f_i,
$$

where $n_v=3$ in triangles and $n_v=4$ in tetrahedra.

`tmesh` applies a small barycentric tolerance (minimum $10^{-8}$) near edges/faces for robustness.

## Self-intersection detection

`find_self_intersections(&surface)` detects intersecting triangles in a 3D boundary
mesh. It returns the list of element index pairs `(i, j)` with `i < j`; an empty list
means that no intersection was detected.

### Algorithm

1. Compute the axis-aligned bounding box (AABB) of every triangle.
2. Build a bounding volume hierarchy (BVH) over these boxes with `parry`.
3. For every triangle, query the BVH for boxes that overlap its own box. This is the
   broad phase: it only keeps candidate pairs.
4. Discard pairs with `i >= j`, so each unordered pair is tested once.
5. Discard pairs that share at least one vertex. Such triangles are topological
   neighbours and touch along an edge or at a vertex by construction.
6. Run the exact triangle-triangle intersection test of `parry` on the remaining
   candidates and keep the intersecting pairs.

The broad phase makes the cost close to linear in the number of triangles for
well-distributed meshes, while the exact test is only applied to the few pairs whose
boxes overlap.

### Limitations

- The input must be a `BoundaryMesh3d`, i.e. a linear triangle mesh. For a volume mesh,
  extract its boundary first with `Mesh::boundary::<BoundaryMesh3d>()`, which returns
  the surface mesh together with the parent face ids.
- Triangles sharing a vertex are never reported. A surface folded onto itself along an
  edge, or two triangles that overlap while sharing a vertex, are therefore missed.
  This choice avoids reporting all the regular adjacencies of the mesh.
- The test is geometric and ignores tags, orientation, and connectivity beyond the
  shared-vertex filter.

```rust
# extern crate tmesh;
use tmesh::Vert3d;
use tmesh::mesh::{BoundaryMesh3d, Mesh, Triangle};
use tmesh::spatialindex::find_self_intersections;

// First triangle in the plane z = 0.
let verts = vec![
	Vert3d::new(0.0, 0.0, 0.0),
	Vert3d::new(1.0, 0.0, 0.0),
	Vert3d::new(0.0, 1.0, 0.0),
	// Second triangle crossing the first one.
	Vert3d::new(0.2, 0.2, -1.0),
	Vert3d::new(0.2, 0.2, 1.0),
	Vert3d::new(0.2, 0.5, 0.3),
];
let elems = vec![Triangle::new(0, 1, 2), Triangle::new(3, 4, 5)];
let etags = vec![1, 1];

let surface = BoundaryMesh3d::from_vecs(verts, elems, etags, vec![], vec![]);
assert_eq!(find_self_intersections(&surface), [(0, 1)]);
```
