# Adjacencies

This page describes which connectivities are stored by default in `tmesh` and which are computed on demand.

## Stored data versus derived adjacency

The core mesh container (`GenericMesh`) stores only primary arrays:

- vertex coordinates,
- element connectivities,
- element tags,
- face connectivities,
- face tags.

Adjacency graphs are computed when needed from these arrays.

## On-demand connectivities

Common derived connectivities include:

- `all_faces()`: canonical face key to adjacent element ids,
- `edges()`: canonical edge key to compact edge id,
- `vertex_to_elems()`: CSR vertex-to-element graph,
- `vertex_to_vertices()`: CSR vertex-to-vertex graph,
- `element_pairs(all_faces)`: element-neighbor pairs across internal faces.

This keeps storage compact and derived connectivities reproducible at any time.

Rust example:

```rust
# extern crate tmesh;
use tmesh::mesh::{rectangle_mesh, Mesh, Mesh2d};

let msh = rectangle_mesh::<Mesh2d>(1.0, 10, 2.0, 15);
let all_faces = msh.all_faces();
let edges = msh.edges();
let v2v = msh.vertex_to_vertices();
let v2e = msh.vertex_to_elems();
let pairs = msh.element_pairs(&all_faces);
let _ = (all_faces, edges, v2v, v2e, pairs);
```

## Dual connectivity

`tmesh` also builds vertex-centered dual structures from simplex primals.

For triangular and tetrahedral meshes, dual construction combines:

- primal edge/face enumeration,
- selected cell-center definitions (median or Barth variants),
- oriented incidence assembly per primal vertex,
- filtering of degenerate contributions.

Outputs can be used either as:

- edge-based normal data for FV-like formulations,
- explicit polygonal/polyhedral dual meshes.

Detailed construction algorithms are documented in [Dual meshes](dual.md).

## Extrusion and unextrusion connectivity transforms

2D triangular meshes can be extruded into 3D prism meshes:

- duplicate vertices at two layers,
- map each triangle to one prism,
- create top/bottom triangular caps and side quadrilateral faces.

The inverse unextrusion extracts the base layer and reconstructs 2D connectivity, followed by orientation fixes.

The algorithms, tag conventions, and orientation rules are detailed in
[Extrusion](extrusion.md).

## See also

- [Tags](tags.md) for face-tag semantics coupled to adjacency inference,
- [Non-manifold meshes](non_manifold.md) for multi-adjacent-face constraints,
- [Dual meshes](dual.md) for dual control-volume construction details,
- [Gradient and Hessian computations](gradients.md) for operators that consume vertex/element adjacency.
