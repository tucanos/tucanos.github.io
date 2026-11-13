# Integration

`tmesh` provides generic integration and $L^2$ norm operators for fields defined at
mesh vertices.

`integrate(f, op)` computes

$$
\int_\Omega \mathrm{op}\left(f_h(x)\right) \,\mathrm{d}x,
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
> **Integration panics**
>
> | Message | Meaning / action |
> | --- | --- |
> | `index out of bounds: the len is {len} but the index is {i}` | `integrate` or `norm` received a field with fewer than `n_verts()` values. |
