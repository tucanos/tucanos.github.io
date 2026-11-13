# Subdivision of canonical elements

When input cells are not simplices, `tmesh` converts them to simplices using deterministic rules based on [*How to Subdivide Pyramids, Prisms and Hexahedra into Tetrahedra.*, Julien Dompierre, Paul Labbe, Marie-Gabrielle Vallet, Ricardo Camarero.](https://www.researchgate.net/publication/221561839_How_to_Subdivide_Pyramids_Prisms_and_Hexahedra_into_Tetrahedra)

Implemented conversions:

- quadrangle $\rightarrow$ 2 triangles,
- pyramid $\rightarrow$ 2 tetrahedra,
- prism $\rightarrow$ 3 tetrahedra,
- hexahedron $\rightarrow$ 5 or 6 tetrahedra.

Internal diagonals are selected from local vertex ordering and edge comparisons so neighboring cells choose compatible splits.

> [!WARNING]
> This subdivision path is implemented for **linear** canonical elements only. High-order (quadratic or above) non-simplex elements are not subdivided by these routines.

> [!WARNING]
> Elements with non-planar quadrangle faces may produce invalid tetrahedra after subdivision.

## Example: create a mesh from one hexahedron and one pyramid

Rust:

```rust
# extern crate tmesh;
use tmesh::{Tag, Vert3d};
use tmesh::mesh::{Hexahedron, Mesh, Mesh3d, Pyramid};

let mut msh = Mesh3d::empty();

// 8 vertices of the unit cube + 1 apex for the pyramid.
let verts = vec![
	Vert3d::new(0.0, 0.0, 0.0), // 0
	Vert3d::new(1.0, 0.0, 0.0), // 1
	Vert3d::new(1.0, 1.0, 0.0), // 2
	Vert3d::new(0.0, 1.0, 0.0), // 3
	Vert3d::new(0.0, 0.0, 1.0), // 4
	Vert3d::new(1.0, 0.0, 1.0), // 5
	Vert3d::new(1.0, 1.0, 1.0), // 6
	Vert3d::new(0.0, 1.0, 1.0), // 7
	Vert3d::new(0.5, 0.5, 1.8), // 8 (pyramid apex)
];
msh.add_verts(verts.into_iter());

// Add one hexahedron (split into 5 or 6 tetrahedra depending on diagonal choice).
let hex = Hexahedron::new([0, 1, 2, 3, 4, 5, 6, 7]);
let _hex_parent_ids = msh.add_hexahedra(std::iter::once(hex), std::iter::once(1 as Tag));

// Add one pyramid (split into exactly 2 tetrahedra).
let pyr = Pyramid::new(4, 5, 6, 7, 8);
msh.add_pyramids(std::iter::once(pyr), std::iter::once(2 as Tag));

assert!(msh.n_elems() >= 7 && msh.n_elems() <= 8); // (5 or 6) + 2
```

Python (`pytucanos`):

```python
import numpy as np
from pytucanos import Mesh3d

msh = Mesh3d.empty()

verts = np.array([
	[0.0, 0.0, 0.0],  # 0
	[1.0, 0.0, 0.0],  # 1
	[1.0, 1.0, 0.0],  # 2
	[0.0, 1.0, 0.0],  # 3
	[0.0, 0.0, 1.0],  # 4
	[1.0, 0.0, 1.0],  # 5
	[1.0, 1.0, 1.0],  # 6
	[0.0, 1.0, 1.0],  # 7
	[0.5, 0.5, 1.8],  # 8 (pyramid apex)
], dtype=float)
msh.add_verts(verts)

hexa = np.array([[0, 1, 2, 3, 4, 5, 6, 7]], dtype=np.int64)
hexa_tags = np.array([1], dtype=np.int16)
_hex_parent_ids = msh.add_hexahedra(hexa, hexa_tags)

pyram = np.array([[4, 5, 6, 7, 8]], dtype=np.int64)
pyram_tags = np.array([2], dtype=np.int16)
msh.add_pyramids(pyram, pyram_tags)

assert 7 <= msh.n_elems() <= 8  # (5 or 6) + 2
```
