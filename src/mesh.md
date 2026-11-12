# Mesh

Mesh handling is implemented primarily in `tmesh` through the `Mesh` trait. It relies on consistent
- vertex coordinates,
- element connectivities,
- face connectivities,
- element tags,
- face tags.

Elements and faces are simplices that can be linear or quadratic. See
[Mesh generation](mesh_generation.md) for the built-in structured and curved
generators and the STL reader, [Mesh I/O](io.md) for reading and writing `.meshb`,
`.solb`, and `.vtu` files, [Extrusion](extrusion.md) for extruding a 2D mesh, and
[Isosurface splitting](isosurface.md) for cutting a mesh along a level set.

## Example: 2D mesh with 3 triangles

Rust:

```rust
# extern crate tmesh;
use tmesh::{Tag, Vert2d};
use tmesh::mesh::{Edge, Mesh, Mesh2d, Triangle};

let verts = vec![
	Vert2d::new(0.0, 0.0),
	Vert2d::new(1.0, 0.0),
	Vert2d::new(2.0, 0.0),
	Vert2d::new(1.0, 1.0),
	Vert2d::new(2.0, 1.0),
];

let elems = vec![
	Triangle::new(0, 1, 3),
	Triangle::new(1, 2, 3),
	Triangle::new(2, 4, 3),
];
let etags: Vec<Tag> = vec![1, 1, 1];

// Boundary edges only (internal edges are not listed here).
let faces = vec![
	Edge::new(0, 1),
	Edge::new(1, 2),
	Edge::new(2, 4),
	Edge::new(4, 3),
	Edge::new(3, 0),
];
let ftags: Vec<Tag> = vec![10, 10, 10, 10, 10];

let msh = Mesh2d::from_vecs(verts, elems, etags, faces, ftags);
assert_eq!(msh.n_elems(), 3);
```

Python (`pytucanos`):

```python
import numpy as np
from pytucanos import Mesh2d

coords = np.array([
	[0.0, 0.0],
	[1.0, 0.0],
	[2.0, 0.0],
	[1.0, 1.0],
	[2.0, 1.0],
], dtype=float)

elems = np.array([
	[0, 1, 3],
	[1, 2, 3],
	[2, 4, 3],
], dtype=np.int64)
etags = np.array([1, 1, 1], dtype=np.int16)

# Boundary edges only.
faces = np.array([
	[0, 1],
	[1, 2],
	[2, 4],
	[4, 3],
	[3, 0],
], dtype=np.int64)
ftags = np.array([10, 10, 10, 10, 10], dtype=np.int16)

msh = Mesh2d(coords, elems, etags, faces, ftags)
assert msh.n_elems() == 3
```

## Partitioning

Partitioning, CSR graphs, and reordering are described in
[Graphs and partitioning](graphs_partitioning.md).
