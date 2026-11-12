# Mesh quality and integration

`tmesh` exposes a small set of geometric quality indicators, plus generic integration
and $L^2$ norm operators for fields defined at the mesh vertices. These indicators are
descriptive: none of them modifies the mesh and none of them returns an error. The only
failure mode is the quadratic-element `gamma` panic listed at the end of this page.

## Face skewness

`face_skewnesses(&all_faces)` iterates over the **internal** faces of the mesh (faces
shared by exactly two elements) and returns `(i0, i1, s)`, where `i0` and `i1` are the
two adjacent elements and $s$ is the skewness of their shared face:

$$
s = \frac{d\left(c_f,\ \mathrm{line}(c_0, c_1)\right)}{|c_1 - c_0|},
$$

with $c_0$ and $c_1$ the centers of the two elements and $c_f$ the center of the shared
face. The denominator is the distance between the two element centers, so $s$ is
dimensionless:

- $s = 0$ when the face center lies on the line joining the element centers, as in an
  orthogonal structured mesh,
- $s$ grows when the face is off-center with respect to the two cells,
- $s$ can exceed $1$ for strongly distorted meshes.

In 2D, a face is an edge and $c_f$ is its midpoint. The `all_faces` argument must be
the face connectivity of the mesh, typically obtained with `Mesh::all_faces()`. Faces
with a number of adjacent elements different from two are skipped.

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 3, 1.0, 3);
let all_faces = msh.all_faces();

let skewnesses: Vec<_> = msh.face_skewnesses(&all_faces).collect();

assert_eq!(skewnesses.len(), 8); // internal edges of the 2 x 2 cell grid
assert!(skewnesses.iter().all(|&(_, _, s)| s.abs() < 1e-6));
```

## Edge length ratio

`edge_length_ratios()` returns, for every element, the ratio between its longest and
its shortest edge:

$$
r_K = \frac{\max_{e \subset K} |e|}{\min_{e \subset K} |e|} \geq 1.
$$

The ratio is $1$ for an equilateral simplex and $\sqrt{2}$ for a right isosceles
triangle. It detects stretched elements, but not all bad elements: a flat tetrahedron
can have almost equal edge lengths.

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 3, 1.0, 3);

// Every triangle of the structured grid is right isosceles.
let sqrt2 = 2.0_f64.sqrt();
assert!(msh
	.edge_length_ratios()
	.all(|r| (r - sqrt2).abs() < 1e-12));
```

## Element quality

`elem_gammas()` returns the geometric quality $\gamma$ of every element, normalized so
that the regular simplex has $\gamma = 1$ and degenerate elements have
$\gamma \to 0$:

| Element | $\gamma$ |
| --- | --- |
| `Edge` | $1$ (always) |
| `Triangle` | $\dfrac{4 \sin A \sin B \sin C}{\sin A + \sin B + \sin C} = \dfrac{2 r_{\mathrm{in}}}{R}$ |
| `Tetrahedron` | $\dfrac{3 r_{\mathrm{in}}}{R}$ |

Here $A, B, C$ are the triangle angles, $r_{\mathrm{in}}$ is the inscribed radius, and
$R$ is the circumradius. For a tetrahedron, $r_{\mathrm{in}} = 3V/S$ where $V$ is the
volume and $S$ the total area of the four faces. The circumradius computation is
ill-conditioned for nearly coplanar vertices; the implementation returns $0$ in that
case (and for zero-volume tetrahedra).

This indicator is useful to locate slivers and needles, especially in 3D, where the
edge ratio alone does not see flat elements. It is only implemented for linear
elements: quadratic triangles and tetrahedra do not provide a `gamma` and trigger an
`unreachable` panic.

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 3, 1.0, 3);

// gamma = 2 (sqrt(2) - 1) for a right isosceles triangle.
let expected = 2.0 * (2.0_f64.sqrt() - 1.0);
assert!(msh
	.elem_gammas()
	.all(|g| (g - expected).abs() < 1e-12));
```

## Integration and norms

`integrate(f, op)` computes

$$
\int_\Omega \mathrm{op}\left(f_h(x)\right) \mathrm{d}x,
$$

where $f$ is a field given at the mesh vertices and $f_h$ is its interpolant on each
element: piecewise linear on a linear mesh and piecewise quadratic on a quadratic mesh.
`norm(f)` is the induced $L^2$ norm:

$$
\|f\|_{L^2(\Omega)} = \left( \int_\Omega f_h(x)^2 \,\mathrm{d}x \right)^{1/2}.
$$

For a boundary mesh, the integrals are taken over the elements of the boundary (edge
lengths in 2D, triangle areas in 3D), not over the enclosed volume.

The element quadratures are:

| Element | Rule | Polynomial exactness |
| --- | --- | --- |
| `Edge` | 3-point Gauss | degree 5 |
| `Triangle` | 3-point (midpoint subdivision) | degree 2 |
| `Tetrahedron` | 4-point | degree 2 |

Because $f_h^2$ is a polynomial of degree 2 on each linear simplex, `norm` is **exact**
for linear meshes. It is not exact for quadratic meshes, where the geometry map is
quadratic and the integrand has a higher degree.

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 10, 1.0, 10);
let f: Vec<f64> = msh.verts().map(|v| v[0]).collect();

// Area of the unit square.
assert!((msh.integrate(&f, |_| 1.0) - 1.0).abs() < 1e-12);
// Integral of x.
assert!((msh.integrate(&f, |x| x) - 0.5).abs() < 1e-12);
// L2 norm of x on the unit square.
assert!((msh.norm(&f) - 1.0 / 3.0_f64.sqrt()).abs() < 1e-12);
```

The field slice is indexed by vertex id; it is not checked. A slice shorter than
`n_verts()` panics with an index-out-of-bounds error, and extra trailing values are
ignored. The field is assumed to be finite; non-finite values propagate to the result.

> [!TIP]
> **Quality and integration panics**
>
> | Message | Meaning / action |
> | --- | --- |
> | `index out of bounds: the len is {len} but the index is {i}` | `integrate` or `norm` received a field with fewer than `n_verts()` values. |
> | `internal error: entered unreachable code` | `elem_gammas` was called on a quadratic mesh; `gamma` is only implemented for linear triangles and tetrahedra. |
>
> `face_skewnesses`, `edge_length_ratios`, and `elem_gammas` return iterators and do not
> allocate a result array; they are usually consumed with `collect`, `map`, `min`, or
> `max`. The face skewness requires the full face connectivity of the mesh, not a
> filtered or submesh connectivity.
