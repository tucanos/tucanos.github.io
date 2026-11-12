# Mesh generation

`tmesh` provides structured generators for canonical geometries and a reader for
surface meshes in STL format. There is no general-purpose Delaunay or advancing-front
mesher: the generators are meant to bootstrap examples, benchmarks, and simple
computational domains, while real geometries are read from `.meshb` or `.stl` files
(see [Mesh I/O](io.md)).

All generators live in `tmesh::mesh` and are generic over the mesh type, so the same
function can produce a linear or a quadratic mesh.

## Structured grids

A structured grid is described by the number of vertices in each direction and by the
vertex coordinates along each axis.

- `nonuniform_rectangle_mesh(&x, &y)` uses the slices `x` and `y` as vertex coordinates.
- `nonuniform_box_mesh(&x, &y, &z)` does the same in 3D.
- `rectangle_mesh(lx, nx, ly, ny)` builds uniform coordinates
  $x_i = i\,l_x/(n_x-1)$ and $y_j = j\,l_y/(n_y-1)$, then calls the nonuniform version.
- `box_mesh(lx, nx, ly, ny, lz, nz)` is the 3D equivalent.

The counts `nx`, `ny`, and `nz` are **numbers of vertices**, not numbers of cells: a
`nx` by `ny` grid contains $(n_x-1)(n_y-1)$ cells. The vertex ordering is
$i + j\,n_x$ in 2D and $i + j\,n_x + k\,n_x n_y$ in 3D.

### 2D grids

Each grid cell is a quadrangle split into two triangles along the diagonal joining its
lowest and highest vertex indices, $(i,j)$ to $(i+1,j+1)$. All element tags are `1`.
Boundary edges are oriented consistently and tagged as follows:

| Face tag | Boundary |
| --- | --- |
| `1` | $y = 0$ |
| `2` | $x = l_x$ |
| `3` | $y = l_y$ |
| `4` | $x = 0$ |

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 4, 2.0, 3);

assert_eq!(msh.n_verts(), 12); // 4 x 3 vertices
assert_eq!(msh.n_elems(), 12); // 3 x 2 cells, 2 triangles each
assert_eq!(msh.n_faces(), 10); // 2 x 4 tag-1/3 edges + 2 x 3 tag-2/4 edges
assert!((msh.vol() - 2.0).abs() < 1e-12);

let mut tags: Vec<_> = msh.ftags().collect();
tags.sort_unstable();
tags.dedup();
assert_eq!(tags, [1, 2, 3, 4]);
```

A nonuniform grid is built directly from the coordinate slices:

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, nonuniform_rectangle_mesh};

let x = [0.0, 0.1, 0.5, 1.0];
let y = [0.0, 1.0];
let msh = nonuniform_rectangle_mesh::<Mesh2d>(&x, &y);

assert_eq!(msh.n_verts(), 8);
assert_eq!(msh.n_elems(), 6);
```

`pytucanos` exposes the same generator with coordinate arrays:

```python
import numpy as np
from pytucanos import Mesh2d

x = np.linspace(0.0, 1.0, 4)
y = np.linspace(0.0, 2.0, 3)
msh = Mesh2d.rectangle_mesh(x, y)

assert msh.n_verts() == 12
assert msh.n_elems() == 12
```

### 3D grids

Each hexahedron is split into five or six tetrahedra, in a way that is consistent
across neighbouring hexahedra: the internal diagonals are selected from the local
vertex ordering and edge comparisons so that neighbouring cells choose compatible
splits (see [Subdivision of canonical elements](subdivision.md)).

All element tags are `1`. Boundary quadrangles are tagged as follows:

| Face tag | Boundary |
| --- | --- |
| `1` | $z = 0$ |
| `2` | $z = l_z$ |
| `3` | $y = 0$ |
| `4` | $y = l_y$ |
| `5` | $x = 0$ |
| `6` | $x = l_x$ |

The 2D grid only fixes the orientation of its boundary edges (`fix_faces_orientation`);
the element orientation is left untouched. Both grid generators therefore assume that
the coordinate slices are increasing. A decreasing slice reverses the orientation of
the incident cells; run `Mesh::fix()` (see [Validity](validity.md)) if the input
coordinates are not monotone.

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh3d, box_mesh};

let msh = box_mesh::<Mesh3d>(1.0, 2, 1.0, 2, 1.0, 2);

assert_eq!(msh.n_verts(), 8); // 2 x 2 x 2 vertices
assert_eq!(msh.n_elems(), 6); // one hexahedron, split into 6 tetrahedra
assert_eq!(msh.n_faces(), 12); // 6 quadrangles, split into 12 triangles
assert!((msh.vol() - 1.0).abs() < 1e-12);
```

## Curved shapes

Curved domains are obtained by splitting a coarse polygonal or polyhedral shape and
projecting vertices onto the exact geometry. The integer argument `n` is a **number of
uniform refinement levels**, not a target edge length: each level splits every element
into $2^d$ elements (2 edges, 4 triangles, or 8 tetrahedra).

| Generator | Starting shape | Projected vertices | `n` |
| --- | --- | --- | --- |
| `circle_mesh(r, n)` | regular $n$-gon | none | number of edges |
| `disk_mesh(n)` | square of side $1/\sqrt{2}$ | boundary vertices onto the circle $r = 1/2$ | levels |
| `sphere_mesh(r, n)` | octahedron surface | all vertices onto the sphere | levels |
| `ball_mesh(r, n)` | octahedron + center | boundary vertices onto the sphere | levels |
| `quadratic_circle_mesh(r, n)` | `circle_mesh` | all nodes | edge count |
| `quadratic_disk_mesh(n)` | `disk_mesh` | boundary nodes | levels |
| `quadratic_sphere_mesh(r, n)` | `sphere_mesh` | all nodes | levels |
| `quadratic_ball_mesh(r, n)` | `ball_mesh` | boundary nodes | levels |

`circle_mesh` and `quadratic_circle_mesh` use `n` as the number of edges of the
polygon; all other generators use it as a refinement count.

> [!WARNING]
> `ball_mesh` and `sphere_mesh` apply the radial projection inside the refinement loop,
> so `n = 0` returns the unit shape regardless of `r`. The quadratic variants always
> project and do not have this limitation.

### Boundary meshes

`circle_mesh` produces a closed 1D boundary mesh of a circle of radius `r` in the plane,
with counter-clockwise edges and tag `1`:

```rust
# extern crate tmesh;
use tmesh::mesh::{BoundaryMesh2d, Mesh, circle_mesh};

let circle = circle_mesh::<BoundaryMesh2d>(1.0, 8);

assert_eq!(circle.n_verts(), 8);
assert_eq!(circle.n_elems(), 8);
// Length of the inscribed octagon.
assert!((circle.vol() - 6.122_934_917_841_436).abs() < 1e-12);
```

`quadratic_circle_mesh` starts from the same polygon, adds one mid-edge node per edge,
and projects every node onto the circle. Its edges are circular arcs, so the length
converges much faster (exactly for a circle):

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, QuadraticBoundaryMesh2d, quadratic_circle_mesh};

let circle = quadratic_circle_mesh::<QuadraticBoundaryMesh2d>(1.0, 8);

assert_eq!(circle.n_verts(), 16); // 8 corners + 8 mid-edge nodes
assert_eq!(circle.n_elems(), 8);
assert!((circle.vol() - 6.280_809_432_296_553).abs() < 1e-12);
```

`sphere_mesh` produces a closed triangulated surface of a sphere in 3D, with tag `1` on
all triangles. Starting from the eight faces of an octahedron, every refinement level
splits each triangle into four and projects all vertices onto the sphere:

```rust
# extern crate tmesh;
use tmesh::mesh::{BoundaryMesh3d, Mesh, sphere_mesh};

let sphere = sphere_mesh::<BoundaryMesh3d>(1.0, 2);

assert_eq!(sphere.n_verts(), 66);
assert_eq!(sphere.n_elems(), 128); // 8 * 4^2
// Area converges to 4 pi.
assert!((sphere.vol() - 11.954_891_630_762_754).abs() < 1e-12);
```

### Volumic meshes

`disk_mesh` starts from the two triangles of a square of side $1/\sqrt{2}$ inscribed in
the circle of radius $1/2$, refines them, and projects the boundary vertices onto that
circle:

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, disk_mesh};

let disk = disk_mesh::<Mesh2d>(2);

assert_eq!(disk.n_verts(), 25); // (2^2 + 1)^2
assert_eq!(disk.n_elems(), 32); // 2 * 4^2
// Area converges to pi / 4.
assert!((disk.vol() - 0.765_366_864_730_179_7).abs() < 1e-12);
```

`ball_mesh` starts from the eight tetrahedra of an octahedron plus its centre, refines
them, and projects the boundary vertices onto the sphere:

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh3d, ball_mesh};

let ball = ball_mesh::<Mesh3d>(1.0, 1);

assert_eq!(ball.n_verts(), 25);
assert_eq!(ball.n_elems(), 64); // 8 * 8
// Volume converges to 4 pi / 3.
assert!((ball.vol() - 2.942_809_041_582_064).abs() < 1e-12);
```

### Quadratic variants

The quadratic generators require the quadratic mesh type as the type parameter. They
first build the linear mesh, add mid-edge nodes, and then project the relevant nodes:

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, QuadraticMesh2d, quadratic_disk_mesh};

let disk = quadratic_disk_mesh::<QuadraticMesh2d>(2);

assert_eq!(disk.n_verts(), 81);
assert_eq!(disk.n_elems(), 32);
// Quadratic triangles represent the circle almost exactly.
assert!((disk.vol() - 0.785_359_429_175_956_9).abs() < 1e-12);
```

## Reading STL surfaces

`read_stl::<M>(file_name)` reads an STL file into a 3D triangle mesh, typically
`BoundaryMesh3d`. Both ASCII and binary STL are handled by the `stl_io` library.

The returned mesh:

- converts the single-precision STL coordinates to `f64`,
- merges vertices whose three coordinates are bit-identical (the `stl_io` reader
  identifies vertices by the bit pattern of their `f32` coordinates, so `0.0` and
  `-0.0` are distinct),
- tags every triangle with `1`,
- has **no explicit face array**.

Because the face array is empty, `Mesh::check()` only passes when the surface is
closed: every edge is then shared by two triangles, so there is no boundary face to
tag. For an open surface, call `Mesh::fix()` to create and orient the boundary edges
before using the mesh with routines that require valid tags.

No orientation fix is applied to the triangles. Run `Mesh::fix()` (or check the
orientation of your file) if the facets may be inconsistently oriented.

```rust
# extern crate tmesh;
use std::io::Write;
use tmesh::mesh::{BoundaryMesh3d, Mesh, read_stl};

let stl = "solid one_triangle\n\
facet normal 0 0 1\n\
outer loop\n\
vertex 0 0 0\n\
vertex 1 0 0\n\
vertex 0 1 0\n\
endloop\n\
endfacet\n\
endsolid one_triangle\n";

let path = std::env::temp_dir().join("one_triangle.stl");
std::fs::File::create(&path)
	.unwrap()
	.write_all(stl.as_bytes())
	.unwrap();

let mut surface = read_stl::<BoundaryMesh3d>(path.to_str().unwrap()).unwrap();
assert_eq!(surface.n_verts(), 3);
assert_eq!(surface.n_elems(), 1);
assert!(surface.etags().all(|t| t == 1));

// Add and orient the boundary edges.
surface.fix().unwrap();
```

## Errors and panics

> [!TIP]
> **Mesh generation error and panic messages**
>
> | Message | Meaning / action |
> | --- | --- |
> | ``called `Result::unwrap()` on an `Err` value: Os { code: 2, kind: NotFound, ... }`` | `read_stl` panics when the file cannot be opened. Check the path. |
> | ``called `Result::unwrap()` on an `Err` value: ...`` (STL parse error) | `read_stl` panics when `stl_io` cannot parse the file. The `Debug` form of the `stl_io` error indicates the failing record; check that the file is a valid ASCII or binary STL. |
>
> `read_stl` never returns an `Err`: both the file opening and the parsing are
> unwrapped, so all I/O and format failures are panics.
>
> The structured generators do not validate their inputs. In particular:
>
> - `rectangle_mesh` and `box_mesh` require at least two vertices per direction. With a
>   count of `0` the `n - 1` capacity computation underflows (panic in debug builds, huge
>   allocation in release builds); with a count of `1` the spacing is infinite and the
>   coordinates become `NaN`.
> - `nonuniform_rectangle_mesh` and `nonuniform_box_mesh` require at least two coordinate
>   values per axis: an empty slice triggers the same underflow. A single value yields a
>   mesh with no element. The values should be increasing, no check is performed, and the
>   element orientation is not fixed.
> - `circle_mesh` requires a positive number of edges; `n = 1` creates a self-loop
>   edge and `n = 0` an empty mesh.
