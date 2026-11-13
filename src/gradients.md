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
	\widehat g_K = \sum_{f \in \partial K} u_{v(f)}\,\mathbf n_f,
$$

where $\mathbf n_f$ is the geometric face normal used in the implementation and
$v(f)$ is the vertex opposite face $f$ (same simplex convention as in [Simplices](simplices.md)).

Then a dimension-dependent normalization is applied when projecting to vertices:

$$
\nabla u_i
=
\frac{\sum_{K\ni i} \widehat g_K}{c_D\,\bar V_i},
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

### Local approximation

The idea is to approximate the variations of a scalar field $u$ on a neighborhood $N(i)$ of vertex $i$ using a local polynomial expansion around $\mathbf{x}_i$. Depending on whether only the gradient or both the gradient and the Hessian are required, either a **linear** or a **quadratic** model is used:

- **Linear least-squares model (gradient reconstruction):**
  $$
  \widehat u_j = \alpha + \sum_{m=1}^{d} a_m \delta_{j,m}
  $$
  with unknown coefficient vector
  $$
  X \equiv \left(\alpha,\, a_1, \dots, a_d\right)^T \in \mathbb{R}^{n_p}, \qquad n_p = 1 + d
  $$
  ($n_p = 3$ in 2D and $n_p = 4$ in 3D). Comparing with the first-order Taylor expansion $u(\mathbf{x}_j) \approx u(\mathbf{x}_i) + \nabla u(\mathbf{x}_i)\cdot\boldsymbol{\delta}_j$, the gradient at $\mathbf{x}_i$ is given directly by
  $$
  \frac{\partial u}{\partial x_m}(\mathbf{x}_i) = a_m, \qquad m \in \{1, \dots, d\}.
  $$

- **Quadratic least-squares model (gradient and Hessian reconstruction):**
  $$
  \widehat u_j = \alpha + \sum_{m=1}^{d} a_m \delta_{j,m} + \sum_{1 \le m \le n \le d} a_{m,n} \delta_{j,m}\delta_{j,n}
  $$
  with unknown coefficient vector
  $$
  X \equiv \left(\alpha,\, a_1, \dots, a_d,\, a_{1,1}, \dots, a_{d,d},\, a_{1,2}, \dots\right)^T \in \mathbb{R}^{n_p}, \qquad n_p = 1 + d + \frac{d(d+1)}{2}
  $$
  ($n_p = 6$ in 2D and $n_p = 10$ in 3D). Comparing with the second-order Taylor expansion $u(\mathbf{x}_j) \approx u(\mathbf{x}_i) + \nabla u(\mathbf{x}_i)\cdot\boldsymbol{\delta}_j + \frac{1}{2}\boldsymbol{\delta}_j^T H(\mathbf{x}_i)\boldsymbol{\delta}_j$, both the gradient and the symmetric Hessian at $\mathbf{x}_i$ are recovered via
  $$
  \frac{\partial u}{\partial x_m}(\mathbf{x}_i) = a_m, \qquad
  H_{mn}(\mathbf{x}_i) = \begin{cases}
  2 a_{m,m} & \text{if } m = n, \\
  a_{m,n} & \text{if } m < n,
  \end{cases}
  $$
  with $H_{nm} = H_{mn}$ by symmetry.

In both models, $\boldsymbol{\delta}_j \equiv \mathbf{x}_j - \mathbf{x}_i$ denotes the displacement vector from vertex $i$ to vertex $j \in N(i)$ with Cartesian components $\delta_{j,m}$, and $X \in \mathbb{R}^{n_p}$ is chosen to minimize the weighted squared residual over $N(i)$:
$$
\mathcal{J}(X) = \sum_{j \in N(i)} w_j^2 \left(\widehat u_j - u_j\right)^2
$$
where $w_j > 0$ are stencil weights—typically inverse-distance weights $w_j = \|\boldsymbol{\delta}_j\|^{-p}$ (e.g., $p = 1$ or $p = 2$) for $j \neq i$, with $w_i = 1$ when vertex $i$ itself is included in $N(i)$.

**Stencil size and second-order neighbors:**
- Determining the $n_p$ polynomial coefficients requires **at least $n_p$ vertices** in $N(i)$ (i.e., $|N(i)| \ge n_p$), and in practice strictly more ($|N(i)| > n_p$) so that the least-squares system is overdetermined rather than interpolating discrete noise.
- For the **linear** model ($n_p = 3$ in 2D, $4$ in 3D), the immediate 1-ring of vertex $i$ (first-order neighbors connected by a mesh edge) almost always provides enough vertices.
- For the **quadratic** model ($n_p = 6$ in 2D, $10$ in 3D), the 1-ring frequently contains fewer than $n_p$ vertices—especially at boundary vertices, corners, or low-valence interior vertices. In such cases, **second-order neighbors** (the 2-ring, i.e., neighbors of first-order neighbors) must be added to $N(i)$ until $|N(i)| \ge n_p$ (and sufficient geometric spread is achieved).

### Matrix form

Ordering the $m = |N(i)| \ge n_p$ stencil vertices as $j = 1, \dots, m$, minimizing $\mathcal{J}(X)$ takes the linear least-squares matrix form
$$
\min_{X \in \mathbb{R}^{n_p}} \left\| W (A X - b) \right\|_2^2 = \min_{X \in \mathbb{R}^{n_p}} \left\| \widehat A X - \widehat b \right\|_2^2
$$
where $\widehat A = W A \in \mathbb{R}^{m \times n_p}$, $\widehat b = W b \in \mathbb{R}^m$, and:
- $W = \operatorname{diag}(w_1, \dots, w_m) \in \mathbb{R}^{m \times m}$ is the diagonal weight matrix,
- $b = (u_1, \dots, u_m)^T \in \mathbb{R}^m$ is the vector of scalar field values on $N(i)$,
- $A \in \mathbb{R}^{m \times n_p}$ is the local Vandermonde design matrix whose $j$-th row evaluates the monomial basis at $\boldsymbol{\delta}_j$.

For the **linear** model, the $j$-th row of $A$ is:
$$
A_{j,:}^{\text{lin, 2D}} = \begin{bmatrix} 1 & \delta_{j,1} & \delta_{j,2} \end{bmatrix}, \qquad
A_{j,:}^{\text{lin, 3D}} = \begin{bmatrix} 1 & \delta_{j,1} & \delta_{j,2} & \delta_{j,3} \end{bmatrix}
$$
For the **quadratic** model, the $j$-th row of $A$ is:
$$
A_{j,:}^{\text{quad, 2D}} = \begin{bmatrix} 1 & \delta_{j,1} & \delta_{j,2} & \delta_{j,1}^2 & \delta_{j,2}^2 & \delta_{j,1}\delta_{j,2} \end{bmatrix}
$$
$$
A_{j,:}^{\text{quad, 3D}} = \begin{bmatrix} 1 & \delta_{j,1} & \delta_{j,2} & \delta_{j,3} & \delta_{j,1}^2 & \delta_{j,2}^2 & \delta_{j,3}^2 & \delta_{j,1}\delta_{j,2} & \delta_{j,1}\delta_{j,3} & \delta_{j,2}\delta_{j,3} \end{bmatrix}
$$

### Solution method & conditioning issues

In anisotropic mesh adaptation, solving $\min_X \|\widehat A X - \widehat b\|_2^2$ directly in physical coordinates suffers from three severe sources of ill-conditioning:
1. **Mesh-size scaling disparity:** Columns of $A$ scale as $O(1)$ (constant $\alpha$), $O(h)$ (linear terms), and $O(h^2)$ (quadratic terms). On refined meshes ($h \ll 1$), column norms differ by orders of magnitude, yielding $\kappa(\widehat A) = O(h^{-1})$ for the linear model and $\kappa(\widehat A) = O(h^{-2})$ for the quadratic model.
2. **Element anisotropy:** On stretched elements with aspect ratio $s = h_{\max}/h_{\min} \gg 1$, displacements along stretched and compressed directions differ by $O(s)$, inflating $\kappa(\widehat A)$ by an additional factor of $O(s)$ (linear) or $O(s^2)$ (quadratic).
3. **Geometric degeneracy:** Near planar boundaries, thin sheets, or collinear neighborhoods, displacements along a normal direction vanish, making columns involving that direction nearly or strictly linearly dependent.

> **Note — Effect of excluding $\alpha$ on conditioning:**
> Fixing $\alpha = u_i$ (i.e., omitting $\alpha$ from the unknowns $X$ and fitting $\widehat u_j - u_i$ over $j \in N(i) \setminus \{i\}$) **significantly reduces the condition number**:
> 1. **Mesh-size scaling:** Removing the $O(1)$ constant column eliminates the largest singular value, lowering the condition number by a full factor of $O(h^{-1})$—from $O(h^{-2})$ down to $O(h^{-1})$ for the quadratic model, and from $O(h^{-1})$ down to $O(1)$ for the linear model.
> 2. **Collinearity with even-degree terms:** Even after normalizing coordinates so that $\|\hat{\boldsymbol{\delta}}_j\| = O(1)$, the constant basis function $1$ and the pure quadratic basis functions $\hat\delta_{j,m}^2$ (which are strictly non-negative and roughly constant on symmetric stencil shells) remain strongly correlated. Omitting $\alpha$ eliminates this near-collinearity between degree-$0$ and degree-$2$ columns.
>
> However, fixing $\alpha = u_i$ forces the reconstructed polynomial to interpolate the discrete nodal value $u_i$ exactly, making the derivatives directly sensitive to point-wise noise or $C^0$ kinks at vertex $i$. Leaving $\alpha$ as a free parameter acts as a local smoothing filter that prevents noise amplification, at the cost of a larger condition number that must be addressed by preconditioning.

#### Normal equations

The classical normal equations set the gradient of $\mathcal{J}(X)$ to zero:
$$
\widehat A^T \widehat A X = \widehat A^T \widehat b \iff A^T W^2 A X = A^T W^2 b
$$
While fast to assemble and solve (a symmetric $n_p \times n_p$ system), forming $\widehat A^T \widehat A$ **squares the condition number**:
$$
\kappa(\widehat A^T \widehat A) = \kappa(\widehat A)^2 = O\!\left(h^{-4} s^4\right) \quad (\text{quadratic})
$$
In IEEE-754 double precision ($\epsilon_{\text{mach}} \approx 2.2 \times 10^{-16}$), all significant digits are lost once $\kappa(\widehat A) \gtrsim 10^8$, making the normal equations unreliable on fine or anisotropic meshes.

#### Thin QR decomposition

Factorizing $\widehat A = Q R$ with $Q \in \mathbb{R}^{m \times n_p}$ ($Q^T Q = I_{n_p}$) and upper-triangular $R \in \mathbb{R}^{n_p \times n_p}$ reduces the least-squares problem to the triangular system
$$
R X = Q^T \widehat b
$$
Thin QR avoids squaring the condition number ($\kappa(R) = \kappa(\widehat A)$). However, without rescaling and rotating the local coordinate frame, $\kappa(\widehat A)$ itself remains $O(h^{-2} s^2)$, and rank-deficient boundary stencils still produce singular or near-singular diagonal entries in $R$.

### Preconditioning via local change of basis

To eliminate both mesh-size scaling ($O(h^{-2})$) and directional stretching ($O(s^2)$), the polynomial is reformulated in a local, dimensionless, isotropic coordinate frame aligned with the principal axes of the stencil $N(i)$.

#### SVD

Collect the weighted relative position vectors into the displacement matrix $D \in \mathbb{R}^{m \times d}$ with rows $D_{j,:} = w_j \boldsymbol{\delta}_j^T$, and compute its thin singular value decomposition:
$$
D = U \Sigma V^T
$$
where:
- $V = [\mathbf{v}_1, \dots, \mathbf{v}_d] \in \mathbb{R}^{d \times d}$ is orthogonal ($V^T V = I_d$) and defines the principal directions of the stencil,
- $\Sigma = \operatorname{diag}(\sigma_1, \dots, \sigma_d)$ with $\sigma_1 \ge \sigma_2 \ge \dots \ge \sigma_d \ge 0$ gives the characteristic weighted stencil length scales along each principal axis.

To prevent division by zero when a stencil is coplanar or collinear ($\sigma_k \ll \sigma_1$), define regularized singular values $\widehat\sigma_k = \max(\sigma_k, \epsilon_{\text{geom}} \sigma_1)$ (with $\epsilon_{\text{geom}} \ll 1$, e.g., $10^{-12}$) and $\widehat\Sigma = \operatorname{diag}(\widehat\sigma_1, \dots, \widehat\sigma_d)$. Defining the linear map $M = V \widehat\Sigma \in \mathbb{R}^{d \times d}$, the normalized local coordinates $\hat{\boldsymbol{\delta}}_j \in \mathbb{R}^d$ are
$$
\hat{\boldsymbol{\delta}}_j = M^{-1} \boldsymbol{\delta}_j = \widehat\Sigma^{-1} V^T \boldsymbol{\delta}_j \iff \boldsymbol{\delta}_j = V \widehat\Sigma \hat{\boldsymbol{\delta}}_j
$$
In this local basis, $\|\hat{\boldsymbol{\delta}}_j\| = O(1)$ isotropically across all active stencil directions.

#### Matrix form

Expressing the local polynomial in the normalized coordinates $\hat{\boldsymbol{\delta}}_j = (\hat\delta_{j,1}, \dots, \hat\delta_{j,d})^T$ (either linear or quadratic) yields the well-scaled local Vandermonde matrix $\hat A \in \mathbb{R}^{m \times n_p}$ (constructed identically to $A$, with $\delta_{j,m}$ replaced by $\hat\delta_{j,m}$) and weighted system
$$
\min_{\hat X \in \mathbb{R}^{n_p}} \left\| \widehat{\hat A} \hat X - \widehat b \right\|_2^2, \qquad \widehat{\hat A} = W \hat A
$$
Because $\boldsymbol{\delta}_j = M \hat{\boldsymbol{\delta}}_j$ is linear, the physical coefficients $X$ and local-basis coefficients $\hat X$ are related by an invertible block-diagonal transformation $T \in \mathbb{R}^{n_p \times n_p}$ such that $\hat A = A T$ and $X = T \hat X$:
- **Constant term:** $\alpha = \hat\alpha$.
- **Gradient (linear and quadratic models):** $\nabla_{\mathbf{x}} u(\mathbf{x}_i) = M^{-T} \hat{\mathbf{a}} = V \widehat\Sigma^{-1} \hat{\mathbf{a}}$.
- **Hessian (quadratic model only):** Forming the symmetric local Hessian $\hat H$ from $\hat X$ ($\hat H_{mm} = 2\hat a_{m,m}$ and $\hat H_{mn} = \hat H_{nm} = \hat a_{m,n}$ for $m < n$), the physical Hessian is
  $$
  H_{\mathbf{x}}(\mathbf{x}_i) = M^{-T} \hat H M^{-1} = V \widehat\Sigma^{-1} \hat H \widehat\Sigma^{-1} V^T
  $$

### SVD preconditioning

Algebraically, transforming the Vandermonde matrix from physical displacements $\boldsymbol{\delta}_j$ to SVD-normalized displacements $\hat{\boldsymbol{\delta}}_j = \widehat\Sigma^{-1} V^T \boldsymbol{\delta}_j$ corresponds to **right-preconditioning** the weighted Vandermonde matrix $\widehat A = W A$ by the block-diagonal polynomial change-of-basis matrix $T$:
$$
\widehat{\hat A} = \widehat A T, \qquad T = \operatorname{diag}\!\left(1,\, M^{-T},\, T_{\text{quad}}\right)
$$
where $T_{\text{quad}}$ maps the local second-order monomials to physical second-order monomials (and is omitted for the linear model).

This SVD-based geometric preconditioning decouples the condition number of the least-squares matrix from both the local mesh size $h$ and the element aspect ratio $s$, reducing $\kappa(\widehat{\hat A})$ to $O(1)$ on non-degenerate stencils.

*(Optional safeguard for rank-deficient stencils)*: When a stencil lies strictly on a lower-dimensional boundary manifold (so $\sigma_d \approx 0$ and all monomial columns involving $\hat\delta_d$ vanish), $\widehat{\hat A}$ has zero columns. In that case, either those degenerate columns are dropped / regularized via a truncated SVD pseudoinverse of $\widehat{\hat A}$ (zeroing singular values below $\tau \bar\sigma_1$), or diagonal regularization is applied so that unobservable normal derivatives default to zero.

### Solving

For each mesh vertex $i$, the complete reconstruction proceeds as follows:

1. **Stencil and weight assembly:** Collect the neighborhood $N(i)$—starting from the 1-ring and adding second-order (or higher-order) neighbors until $|N(i)| \ge n_p$—compute displacements $\boldsymbol{\delta}_j = \mathbf{x}_j - \mathbf{x}_i$, and evaluate the diagonal weights $W = \operatorname{diag}(w_1, \dots, w_m)$.
2. **Local geometric change of basis (SVD preconditioning):** Form the weighted displacement matrix $D \in \mathbb{R}^{m \times d}$, compute its thin SVD $D = U \Sigma V^T$, regularize small singular values to form $\widehat\Sigma$, and map displacements to normalized local coordinates $\hat{\boldsymbol{\delta}}_j = \widehat\Sigma^{-1} V^T \boldsymbol{\delta}_j$.
3. **Preconditioned system assembly and Thin QR solve:** Assemble $\widehat{\hat A} = W \hat A$ and $\widehat b = W b$, compute the thin QR factorization $\widehat{\hat A} = Q R$, and solve the upper-triangular system
   $$
   R \hat X = Q^T \widehat b
   $$
   for the local-basis coefficients $\hat X$ (falling back to a truncated SVD pseudoinverse if $R$ is numerically rank-deficient).
4. **Back-transformation to physical space:** Extract the local gradient $\hat{\mathbf{a}}$ (and, for the quadratic model, the symmetric local Hessian $\hat H$) from $\hat X$, and transform back to global Cartesian coordinates:
   $$
   \nabla_{\mathbf{x}} u(\mathbf{x}_i) = V \widehat\Sigma^{-1} \hat{\mathbf{a}}, \qquad
   H_{\mathbf{x}}(\mathbf{x}_i) = V \widehat\Sigma^{-1} \hat H \widehat\Sigma^{-1} V^T
   $$
   
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
