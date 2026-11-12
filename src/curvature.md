# Curvature

## Curvature tensor: geometric meaning

For a smooth surface, curvature is the differential variation of the normal field inside the tangent plane. This is the (symmetric) shape operator.

In an orthonormal tangent basis $(\mathbf{e}_u, \mathbf{e}_v)$, it is represented by a $2\times2$ symmetric matrix

$$
S = \begin{bmatrix}
k_{uu} & k_{uv} \\
k_{uv} & k_{vv}
\end{bmatrix}.
$$

Its eigenvalues are the principal curvatures $(k_1, k_2)$ and its eigenvectors are the principal directions.

In `tucanos`, the stored output is not the full matrix per element. Instead, the algorithm returns principal direction vectors already scaled by principal curvatures:

- $\mathbf{u} = k_1\,\hat{\mathbf{u}}$,
- $\mathbf{v} = k_2\,\hat{\mathbf{v}}$ (in 3D surface case).

So vector directions encode principal directions, and vector norms encode principal curvatures.

## What is computed in each case

- **2D edge mesh**: one curvature direction per element (`u`).
- **3D triangle mesh**: two principal directions per element (`u`, `v`).


## Discrete algorithm used in tucanos

For quadratic meshes, curvature is obtained from quadratic geometric elements (`QuadraticGEdge` / `QuadraticGTriangle`) directly. For linear meshed, the algorithm described [here](https://gfx.cs.princeton.edu/pubs/Rusinkiewicz_2004_ECA/curvpaper.pdf) is outlined below 

### 1. Vertex normal estimation

Normals are first accumulated from incident element normals and normalized.

Available weights:

- `Volume` (default in production path),
- `VolumeOverEdgeLengthSquared` (mainly used in tests for triangles).

Important assumption: normals should not be aggregated across non-smooth feature boundaries.

### 2. Element-wise tensor fit 

For each element, normal differences across edges are used to estimate curvature coefficients.

- In 2D, this gives a scalar curvature along local tangent direction.
- In 3D triangles, a local tangent frame $(\mathbf{u},\mathbf{v})$ is built and a least-squares system is assembled for $(k_{uu},k_{uv},k_{vv})$.

`tucanos` solves this small system with QR factorization.

### 3. Principal decomposition (3D)

The symmetric matrix

$$
\begin{bmatrix}
k_{uu} & k_{uv} \\
k_{uv} & k_{vv}
\end{bmatrix}
$$

is diagonalized (symmetric eigendecomposition). Eigenvectors define principal directions; eigenvalues define principal curvatures.

### 4. Magnitude bounding

Curvature magnitudes are bounded using the global scale `H_MAX` (through `bound_curvature`) so small/near-zero values are stabilized and extreme values are clipped consistently with remeshing constraints.

### 5. Boundary fix-up

Curvature near patch boundaries is less reliable. `fix_curvature` propagates orientation/magnitude from valid interior neighbors with the same element tag to boundary-adjacent elements.

This yields more consistent tensors around patch borders.

## Notes on robustness and limitations

- Estimation is less accurate near boundaries before fix-up.
- Curvature across non-smooth patch junctions is not meaningful unless features are handled/tagged consistently.
- CAD/meshed geometry consistency and topology tagging quality directly impact curvature quality.
