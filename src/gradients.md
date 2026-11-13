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

The idea is to approximate the variations of a scalar field $u$ on a neighborhood $N(i)$ of vertex $i$ using a local polynomial expansion around $\mathbf{x}_i$. To avoid baseline offsets, field values are shifted by the central vertex value $u_i$, so we fit the differences $\Delta u_j = u_j - u_i$.

Depending on whether only the gradient or both the gradient and the Hessian are required, either a **linear** (`order = 1`) or a **quadratic** (`order = 2`) weighted least-squares model is used:

* **Linear least-squares model (gradient reconstruction):**

$$
\widehat u_j = \alpha + \sum_{m=1}^{d} a_m \delta_{j,m}
$$

with unknown coefficient vector

$$
X \equiv \left(\alpha,\, a_1, \dots, a_d\right)^T \in \mathbb{R}^{n_p}, \qquad n_p = 1 + d
$$

($n_p = 3$ in 2D and $n_p = 4$ in 3D). The physical gradient at $\mathbf{x}_i$ is directly given by $\frac{\partial u}{\partial x_m}(\mathbf{x}_i) = a_m$.

* **Quadratic least-squares model (gradient and Hessian reconstruction):**

$$
\widehat u_j = \alpha + \sum_{m=1}^{d} a_m \delta_{j,m} + \sum_{m=1}^{d} \frac{1}{2} a_{m,m} \delta_{j,m}^2 + \sum_{1 \le m < n \le d} a_{m,n} \delta_{j,m}\delta_{j,n}
$$

with unknown coefficient vector

$$
X \equiv \left(\alpha,\, a_1, \dots, a_d,\, a_{1,1}, \dots, a_{d,d},\, a_{1,2}, \dots\right)^T \in \mathbb{R}^{n_p}, \qquad n_p = 1 + d + \frac{d(d+1)}{2}
$$

($n_p = 6$ in 2D and $n_p = 10$ in 3D). Including the explicit $\frac{1}{2}$ factor on the pure quadratic terms allows both the gradient and the symmetric Hessian at $\mathbf{x}_i$ to be recovered directly without rescaling diagonal entries:

$$
\frac{\partial u}{\partial x_m}(\mathbf{x}_i) = a_m, \qquad H_{mn}(\mathbf{x}_i) = a_{m,n}
$$

with $H_{nm} = H_{mn}$ by symmetry.

**Stencil size and second-order neighbors:**
To determine all $n_p$ polynomial coefficients, the stencil must contain at least $n_p$ neighbors ($m = \vert{}N(i)\vert{} \ge n_p$):

* For the **linear** model, the immediate 1-ring of vertex $i$ is used.

* For the **quadratic** model, the 1-ring frequently contains fewer than $n_p$ vertices (or produces a poorly poised stencil), so **second-order neighbors** (the 2-ring, excluding vertex $i$ itself) are automatically included to ensure an overdetermined system.

In both models, $\boldsymbol{\delta}_j \equiv \mathbf{x}_j - \mathbf{x}_i$ denotes the displacement vector from vertex $i$ to neighbor $j$. The coefficient vector $X$ is chosen to minimize the weighted squared residual over $N(i)$, augmented with a Tikhonov regularization term on the constant correction $\alpha$:

$$
\mathcal{J}(X) = \sum_{j \in N(i)} w_j^2 \left(\widehat u_j - \Delta u_j\right)^2 + \lambda \alpha^2
$$

where:

* $w_j = \Vert{}\boldsymbol{\delta}_j\Vert{}^{-p}$ (`p = weight_exp`, e.g., $p \in \{0, 1, 2\}$) are inverse-distance weights normalized by $w_{\max} = \max_{j \in N(i)} w_j$ so that $\max_j w_j = 1$.

* $\lambda = 2 w_{\max}^2 = 2$ penalizes deviations of $\alpha$ from zero, acting as a central-vertex anchor that prevents ill-posed constant shifts while still allowing $\alpha$ to serve as a local smoothing correction ($u_i^{\text{smooth}} = u_i + \alpha$).

### Matrix form

Ordering the $m = \vert{}N(i)\vert{}$ stencil neighbors as $j = 1, \dots, m$, minimizing $\mathcal{J}(X)$ takes the linear least-squares matrix form

$$
\min_{X \in \mathbb{R}^{n_p}} \left\Vert{} W (A X - \Delta U) \right\Vert{}_2^2 + \lambda \alpha^2 = \min_{X \in \mathbb{R}^{n_p}} \left\Vert{} \widehat A X - \widehat b \right\Vert{}_2^2 + \lambda \alpha^2
$$

where $\widehat A = W A \in \mathbb{R}^{m \times n_p}$, $\widehat b = W \Delta U \in \mathbb{R}^m$, and:

* $W = \operatorname{diag}(w_1, \dots, w_m) \in \mathbb{R}^{m \times m}$ is the diagonal matrix of normalized weights,

* $\Delta U = (\Delta u_1, \dots, \Delta u_m)^T \in \mathbb{R}^m$ is the vector of field differences,

* $A \in \mathbb{R}^{m \times n_p}$ is the local Vandermonde design matrix whose $j$-th row evaluates the polynomial basis $\boldsymbol{\phi}(\boldsymbol{\delta}_j)^T$ at $\boldsymbol{\delta}_j$.

For the quadratic model (`order = 2`), the $j$-th row of $A$ in 2D ($d=2$, $n_p=6$) and 3D ($d=3$, $n_p=10$) is:

$$
A_{j,:}^{\text{2D}} = \begin{bmatrix} 1 & \delta_{j,1} & \delta_{j,2} & \frac{1}{2}\delta_{j,1}^2 & \frac{1}{2}\delta_{j,2}^2 & \delta_{j,1}\delta_{j,2} \end{bmatrix}
$$

$$
A_{j,:}^{\text{3D}} = \begin{bmatrix} 1 & \delta_{j,1} & \delta_{j,2} & \delta_{j,3} & \frac{1}{2}\delta_{j,1}^2 & \frac{1}{2}\delta_{j,2}^2 & \frac{1}{2}\delta_{j,3}^2 & \delta_{j,1}\delta_{j,2} & \delta_{j,2}\delta_{j,3} & \delta_{j,1}\delta_{j,3} \end{bmatrix}
$$

(For `order = 1`, only the first $1 + d$ columns are retained.)

### Solution method & conditioning issues

Solving the weighted least-squares problem directly in physical coordinates suffers from three primary sources of ill-conditioning on adapted meshes:

1. **Mesh-size scaling disparity:** The columns of $A$ scale as $O(1)$ (constant term), $O(h)$ (linear terms), and $O(h^2)$ (quadratic terms). On refined meshes ($h \ll 1$), column norms span multiple orders of magnitude, giving $\kappa(\widehat A) = O(h^{-2})$.

2. **Element anisotropy:** On stretched elements with aspect ratio $s = h_{\max}/h_{\min} \gg 1$ (common in boundary-layer or shock-aligned meshes, where $s$ can reach $10^6$), displacements along stretched and compressed directions differ by $O(s)$, inflating $\kappa(\widehat A)$ by an additional factor of $O(s^2)$.

3. **Geometric degeneracy:** When a stencil is nearly coplanar or collinear (e.g., at boundary vertices or collapsed neighborhoods), displacements along a normal direction vanish, making columns associated with that direction nearly linearly dependent.

> [!NOTE]
> **Excluding $\alpha$:** Setting $\alpha = 0$ (enforcing exact interpolation $\widehat u_i = u_i$ at vertex $i$ and dropping the first column of $A$) removes the $O(1)$ constant column. This reduces the mesh-size condition number by a full factor of $O(h^{-1})$—from $O(h^{-2})$ to $O(h^{-1})$ for the quadratic model, and from $O(h^{-1})$ to $O(1)$ for the linear model—and avoids coupling between the constant column $1$ and the non-negative quadratic columns $\delta_{j,m}^2$. In `tucanos`, $\alpha$ is retained so the same operator can filter noise (`smooth`), while the damping term $\lambda \alpha^2$ ($\lambda = 2 w_{\max}^2$) stabilizes the first column.

#### Normal equations

Setting the gradient of $\mathcal{J}(X)$ with respect to $X$ to zero yields the regularized normal equations:

$$
\left(A^T W^2 A + \Lambda\right) X = A^T W^2 \Delta U, \qquad \Lambda = \operatorname{diag}(\lambda, 0, \dots, 0)
$$

* **Advantage:** The system matrix $N = A^T W^2 A + \Lambda$ is symmetric positive-definite of tiny fixed size $n_p \times n_p$ ($3\times 3$ or $6\times 6$ in 2D; $4\times 4$ or $10\times 10$ in 3D), making it very fast to assemble and factorize via Cholesky decomposition. Moreover, $N^{-1} A^T W^2$ depends only on mesh geometry and can be reused across multiple scalar fields.

* **Conditioning drawback:** Forming $A^T W^2 A$ **squares the condition number**:

$$
\kappa(N) \approx \kappa(\widehat A)^2 = O\!\left(h^{-4} s^4\right)
$$

In IEEE-754 double precision ($\epsilon_{\text{mach}} \approx 2.2 \times 10^{-16}$), all accuracy is lost once $\kappa(\widehat A) \gtrsim 10^8$, causing Cholesky factorization to fail or return corrupted derivatives if applied directly in unscaled physical coordinates.

#### Thin QR decomposition

Factorizing the augmented weighted design matrix $\widehat A_{\text{aug}} = Q R$ with $Q \in \mathbb{R}^{(m+1) \times n_p}$ ($Q^T Q = I_{n_p}$) and upper-triangular $R \in \mathbb{R}^{n_p \times n_p}$ reduces the problem to the triangular system

$$
R X = Q^T \widehat b_{\text{aug}}
$$

* **Advantage:** Thin QR avoids squaring the condition number ($\kappa(R) = \kappa(\widehat A_{\text{aug}})$), making it substantially more accurate than unpreconditioned normal equations.

* **Drawback:** Computing a QR factorization of an $m \times n_p$ matrix at every mesh vertex is more expensive than assembling and Cholesky-factoring an $n_p \times n_p$ normal matrix. Furthermore, without changing to a local isotropic basis, $\kappa(R) = O(h^{-2} s^2)$ still degrades on extreme aspect ratios ($s \sim 10^6$).

### Preconditioning via local change of basis

To combine the speed of the normal equations with numerical stability on arbitrary mesh sizes and extreme anisotropies, the polynomial is reformulated in a local, dimensionless, isotropic coordinate frame aligned with the principal axes of the weighted stencil.

#### SVD

Collect the weighted relative displacement vectors into the $m \times d$ coordinate matrix $X_w \in \mathbb{R}^{m \times d}$ with rows

$$
[X_w]_{j,:} = w_j \boldsymbol{\delta}_j^T
$$

and compute its thin singular value decomposition:

$$
X_w = U \Sigma V^T
$$

where:

* $V = [\mathbf{v}_1, \dots, \mathbf{v}_d] \in \mathbb{R}^{d \times d}$ is orthogonal ($V^T V = I_d$) and defines the principal geometric directions of the stencil,

* $\Sigma = \operatorname{diag}(\sigma_1, \dots, \sigma_d)$ with $\sigma_1 \ge \sigma_2 \ge \dots \ge \sigma_d \ge 0$ captures the characteristic weighted stencil length scales along each principal axis.

To prevent division by zero on flat or degenerate stencils ($\sigma_k \ll \sigma_1$), the inverse singular values are clamped using a relative threshold of $10^{-14}$:

$$
\sigma_k^{-1} = \frac{1}{\max\!\left(\sigma_k,\, \max(\sigma_1, 10^{-20}) \cdot 10^{-14}\right)}
$$

This defines the forward whitening transformation $P \in \mathbb{R}^{d \times d}$ and its inverse transpose $T_{\text{inv}} \in \mathbb{R}^{d \times d}$:

$$
P = \Sigma^{-1} V^T, \qquad T_{\text{inv}} = P^{-T} = V \Sigma^{-1}
$$

Every displacement vector is mapped into the whitened local frame prior to polynomial evaluation:

$$
\boldsymbol{\delta}'_j = P \boldsymbol{\delta}_j = \Sigma^{-1} V^T \boldsymbol{\delta}_j
$$

#### Matrix form

Evaluating the unweighted polynomial basis $\boldsymbol{\phi}(\boldsymbol{\delta}'_j)$ at the transformed coordinates $\boldsymbol{\delta}'_j$ yields the preconditioned Vandermonde matrix $A' \in \mathbb{R}^{m \times n_p}$ and transformed unknown vector $X'$:

$$
\min_{X' \in \mathbb{R}^{n_p}} \left\Vert{} W (A' X' - \Delta U) \right\Vert{}_2^2 + \lambda \alpha^2
$$

Because $\boldsymbol{\delta}'_j = P \boldsymbol{\delta}_j$ is a linear change of coordinates, the local derivatives $\mathbf{g}' = (a'_1, \dots, a'_d)^T$ and symmetric local Hessian $H'$ ($H'_{mn} = a'_{m,n}$) map back to physical Cartesian space via the chain rule:

$$
\nabla_{\mathbf{x}} u(\mathbf{x}_i) = T_{\text{inv}} \mathbf{g}', \qquad H_{\mathbf{x}}(\mathbf{x}_i) = T_{\text{inv}} H' T_{\text{inv}}^T
$$

### SVD preconditioning

Applying the $d \times d$ SVD whitening transform $P = \Sigma^{-1} V^T$ to the coordinates acts as an exact geometric right-preconditioner on the Vandermonde system:

1. **Elimination of** $O(h^{-2})$ **and** $O(s^2)$ **ill-conditioning:** In the whitened frame, the weighted displacements $w_j \boldsymbol{\delta}'_j$ have unit singular values along all non-degenerate principal axes. Consequently, the columns of $W A'$ are $O(1)$ regardless of mesh refinement $h$ or element stretching $s$, keeping $\kappa(W A') = O(1)$ even for aspect ratios as severe as $s = 10^6$.

2. **Safe use of Normal Equations + Cholesky:** Because $\kappa(W A') = O(1)$, squaring the condition number in the preconditioned normal matrix

$$
N = (A')^T W^2 A' + \Lambda
$$

yields $\kappa(N) = \kappa(W A')^2 = O(1)$, well within double-precision limits. This eliminates the need for a more expensive Thin QR factorization or a second SVD on the $m \times n_p$ Vandermonde matrix.

3. **Fallback for boundary or degenerate vertices:** Boundary vertices (`boundary_flag`) and any interior vertices where Cholesky factorization fails due to extreme geometric degeneracy are flagged and subsequently filled by iteratively averaging valid values from their 1-ring neighbors (`fix_not_computed`).

### Solving

For each non-boundary mesh vertex $i$, the reconstruction proceeds as follows:

1. **Neighbor collection:** Gather the 1-ring of vertex $i$ (and expand to the 2-ring, excluding $i$, when `order = 2`), then compute physical displacements $\boldsymbol{\delta}_j = \mathbf{x}_j - \mathbf{x}_i$ and field differences $\Delta u_j = u_j - u_i$.

2. **Weight normalization:** Compute $w_j = \Vert{}\boldsymbol{\delta}_j\Vert{}^{-p}$ and normalize all weights by $w_{\max}$ so that $\max_j w_j = 1$.

3. **Coordinate whitening via SVD:** Assemble the $m \times d$ matrix $X_w$ with rows $w_j \boldsymbol{\delta}_j^T$, compute its thin SVD $X_w = U \Sigma V^T$, clamp small singular values, and form $P = \Sigma^{-1} V^T$ and $T_{\text{inv}} = V \Sigma^{-1}$.

4. **Normal equations assembly:** Map each offset to $\boldsymbol{\delta}'_j = P \boldsymbol{\delta}_j$, evaluate the polynomial basis $\boldsymbol{\phi}(\boldsymbol{\delta}'_j)$, and accumulate the symmetric $n_p \times n_p$ normal matrix and $n_p \times m$ projection matrix:

$$
N = \Lambda + \sum_{j=1}^{m} w_j^2 \boldsymbol{\phi}(\boldsymbol{\delta}'_j) \boldsymbol{\phi}(\boldsymbol{\delta}'_j)^T, \qquad B = \begin{bmatrix} w_1^2 \boldsymbol{\phi}(\boldsymbol{\delta}'_1) & \cdots & w_m^2 \boldsymbol{\phi}(\boldsymbol{\delta}'_m) \end{bmatrix}
$$

with $\Lambda = \operatorname{diag}(2 w_{\max}^2, 0, \dots, 0)$.

5. **Cholesky solve:** Factorize $N = L L^T$ via Cholesky decomposition and solve $N X' = B \Delta U$ for $X'$.

6. **Physical derivative recovery and boundary propagation:** Extract the smoothed correction $\alpha = X'_0$, map $\mathbf{g}'$ and $H'$ back to physical coordinates via $\nabla_{\mathbf{x}} u = T_{\text{inv}} \mathbf{g}'$ and $H_{\mathbf{x}} = T_{\text{inv}} H' T_{\text{inv}}^T$, and propagate valid values to boundary or failed vertices via iterative 1-ring neighbor averaging.


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
