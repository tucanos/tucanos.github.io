# Gradient and Hessian computations

`tmesh` exposes three derivative strategies:

- `LinearLeastSquares(p)`: weighted least-squares gradient (linear model),
- `QuadraticLeastSquares(p)`: weighted least-squares gradient/Hessian (quadratic model),
- `L2Projection`: element-based gradient projection.

## L2 projection gradient

For scalar values $u_i$ at vertices:

$$
\nabla u_i = \frac{\sum_{K \ni i} |K|\,\nabla u_K}{\sum_{K \ni i} |K|}.
$$

Element gradient in L2 projection:

For one simplex element $K$ (triangle in 2D or tetrahedron in 3D), `tmesh` builds
an element-wise gradient from face contributions:

$$
	\tilde g_K = \sum_{f \in \partial K} u_{v(f)}\,\mathbf n_f,
$$

where $\mathbf n_f$ is the geometric face normal used in the implementation and
$v(f)$ is the vertex opposite face $f$ (same simplex convention as in [Simplices](simplices.md)).

Then a dimension-dependent normalization is applied when projecting to vertices:

$$
\nabla u_i
=
\frac{\sum_{K\ni i} \tilde g_K}{c_D\,\bar V_i},
\qquad
\bar V_i = \frac{1}{n_v}\sum_{K\ni i}|K|,
$$

with $n_v$ the number of element vertices ($n_v=3$ in triangles, $n_v=4$ in tetrahedra),
and

$$
c_2=-6,\qquad c_3=-12.
$$

This is equivalent to the volume-weighted expression above, written in the exact scaling form used by `tmesh`.

Implementation details:

- element gradients are assembled from face normals and vertex values,
- accumulation uses vertex-to-element adjacency,
- scaling follows simplex dimension-dependent constants.

For Hessians, L2 projection applies gradient projection twice and symmetrizes second-order tensors; this is a rough approximation for Hessians.

## Weighted least-squares model

Around vertex $\mathbf{x}_i$, for neighbors $j\in N(i)$ define:

$$
\delta^{(j)} = \mathbf{x}_j - \mathbf{x}_i,
\qquad
\Delta f_j = f_j - f_i.
$$

Quadratic local model:

$$
f(\mathbf{x}) \approx \alpha + a_m\delta_m + a_{m,n}\delta_m\delta_n.
$$

First-order mode omits $a_{m,n}$ terms.

Objective minimized at each vertex:

$$
W_0\alpha^2 + \sum_{j\in N(i)} W_j\left(\Delta f_j - \left(\alpha + a_k\delta_k^{(j)} + a_{k,\ell}\delta_k^{(j)}\delta_\ell^{(j)}\right)\right)^2,
$$

with

$$
W_j = \frac{1}{\|\delta^{(j)}\|^p}.
$$

In second-order mode, the basis uses:

- $0.5\,\delta_m^2$ for diagonal quadratic terms,
- $\delta_m\delta_n$ for cross terms.

## Solver alternatives and selected strategy

Common LS backend choices are:

- normal equations + Cholesky,
- thin QR,
- preconditioned normal equations + Cholesky.

At the pinned reference baseline, `tmesh` uses preconditioned normal equations + Cholesky with SVD coordinate preconditioning.

## SVD preconditioning and solve

Let $X_w\in\mathbb{R}^{n\times d}$ have rows $w_j\,\delta^{(j)\top}$ and SVD:

$$
X_w = U\,\Sigma\,V^\top.
$$

Using clipped inverse singular values $\Sigma^{-1}_{\varepsilon}$:

$$
P = \Sigma^{-1}_{\varepsilon}V^\top,
\qquad
T^{-1}=V\,\Sigma^{-1}_{\varepsilon}.
$$

The polynomial basis is evaluated in transformed coordinates
$\tilde\delta^{(j)} = P\,(w_j\delta^{(j)})$.

With transformed design matrix $A$ and diagonal weights $W$:

$$
N = A^\top A,
\qquad
B = A^\top W,
\qquad
N\,x = B\,\Delta f,
$$

and $N\,x=\cdot$ is solved by Cholesky.

## Gradient and Hessian recovery

From solved coefficients $x$:

- extract transformed linear coefficients $g'$,
- assemble transformed symmetric Hessian matrix $H'$ from quadratic slots,
- recover physical derivatives by

$$
\nabla f = T^{-1}g',
\qquad
H = T^{-1}H'\left(T^{-1}\right)^\top.
$$

Hessian component output ordering:

- 2D: $[H_{xx}, H_{yy}, H_{xy}]$,
- 3D: $[H_{xx}, H_{yy}, H_{zz}, H_{xy}, H_{yz}, H_{xz}]$.

## Conditioning and failure recovery

Conditioning controls include:

- distance-based row scaling via $W_j$,
- anchor scaling $W_0 = \sqrt{2}\,\max_j W_j$,
- SVD singular-value clipping,
- conditioning checks before accepting local solves.

If a local solve fails (too few neighbors or poor conditioning), `tmesh` marks the vertex unresolved and recovers values by iterative neighbor averaging.

## Errors and method restrictions

The least-squares backend reports local-solve failures, and the mesh-level
`gradient`/`hessian` entry points temporarily mark those vertices as unresolved
before running the recovery pass.

> [!TIP]
> **Least-squares error messages**
>
> | Message | Meaning / action |
> | --- | --- |
> | `Failed Cholesky decomposition of normal equations` | The preconditioned normal matrix is not positive definite, typically because the neighborhood is degenerate or too few independent neighbors are available. The vertex is marked unresolved and recovered by averaging. |
> | `NaN in gradient computation` | The recovered field still contains `NaN` after the recovery pass. Check the input field for `NaN`/`Inf` and avoid coincident or degenerate vertices. |
> | `Cannot compute the value everywhere` | Too many vertices could not be resolved (for example, boundary-only neighborhoods) and neighbor averaging could not fill them. Use more neighbors or a first-order model. |
>
> The mesh-level `Mesh::gradient` and `Mesh::hessian` methods call `.unwrap()` on
> these results, so in practice the messages appear in a panic rather than in a
> returned `Result`.

> [!TIP]
> **Method restrictions**
>
> Two combinations are rejected with a panic instead of an error:
>
> - `Cannot use LinearLeastSquares to compute the hessian` — the linear model has no
>   second-order terms; use `QuadraticLeastSquares` or `L2Projection`.
> - `Cannot use L2Proj for smoothing` — L2 projection is not a smoothing operator;
>   use one of the least-squares methods for `smooth`.

Rust example:

```rust
# extern crate tmesh;
use tmesh::mesh::{GradientMethod, Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 20, 1.0, 20);
let f: Vec<f64> = msh.verts().map(|x| x[0] * x[0] + x[1]).collect();

let g = msh.gradient(GradientMethod::QuadraticLeastSquares(2), &f);
let h = msh.hessian(GradientMethod::QuadraticLeastSquares(2), &f);
let _ = (g, h);
```

Python example:

```python
from pytucanos.mesh import get_square
from pytucanos import Mesh2d

coords, elems, etags, faces, ftags = get_square()
msh = Mesh2d(coords, elems, etags, faces, ftags).split().split()
x, y = msh.get_verts().T
f = (x * x + y).reshape((-1, 1))

g = msh.gradient(f, weight_exp=2, order=1)
h = msh.hessian(f, weight_exp=2)
```
