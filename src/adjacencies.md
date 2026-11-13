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
