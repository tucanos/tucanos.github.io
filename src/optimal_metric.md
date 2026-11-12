# Optimal metric definition
## Feature based adaptation
Feature-based metric adaptation consists in applying the **multiscale error estimator** in $L^p$ norm to a scalar function in the domain. This function is denoted as a sensor (or feature), and the mesh is adapted to reduce global interpolation error on this feature while capturing multiple scales. In CFD computations, the Mach number is often used because it captures boundary layers, shocks, and wakes.

Let $\mathcal{T}_h$ be a simplex mesh composed of simplex elements $K$. Let $E_{L^p}(\mathcal{T}_h)$ be an error model in $L^p$-norm. The optimal mesh $\mathcal{T}_h^{opt}$ is given by the following minimization problem with constraint $\mathcal{C}(\mathcal{T}_h)=N$:

$$
\mathcal{T}_h^{opt}=\operatorname*{argmin}_{\mathcal{C}(\mathcal{T}_h)=N} {E_{L^p}(\mathcal{T}_h)}
$$

In the continuous mesh theoretical framework, the goal is to find the optimal continuous mesh $\mathcal{M}_{opt}$ which minimizes the given continuous error model $\mathcal{E}_{L^p}$ for a fixed continuous mesh complexity $C(\mathcal{M}) = N$:

$$
\mathcal{M}_{opt}=\operatorname*{argmin}_{\mathcal{C}(\mathcal{M})=N} {\mathcal{E}_{L^p}(\mathcal{M})}
$$

## Multiscale error estimator in $L^p$-norm and feature-based mesh adaptation

Let $w=f(\mathbf{u})$ be a twice-differentiable smooth function of state $\mathbf{u}$. Let $\mathcal{T}_h$ be a simplex mesh on which $w$ is approximated piecewise linearly on elements $K$. The approximation is defined by linear interpolant $\Pi_h w$ on mesh $\mathcal{T}_h$. Then the global interpolation error of $w$ in $\mathcal{T}_h$ is the $L^p$ norm of local interpolation errors:

$$
E_{L^p}(\mathcal{T}_h) = \left(\int_{\Omega}|w-\Pi_hw|^p d\Omega\right)^{\frac{1}{p}}.
$$

In Loseille (2011), existence and uniqueness of continuous linear interpolant $\pi_{\mathcal{M}}w$ are established for continuous mesh $(\mathcal{M}(\mathbf{x}))_{\mathbf{x}\in\Omega}$. From Corollary 3.4, the following continuous interpolation estimate holds in 2D and 3D:

$$
\begin{aligned}
    \forall \mathbf{a}\in\Omega, \quad |w-\pi_{\mathcal{M}}w|(\mathbf{a})&=\frac{1}{8}\text{tr}(\mathcal{M}(\mathbf{a})^{-\frac{1}{2}} |H(\mathbf{a})| \mathcal{M}(\mathbf{a})^{-\frac{1}{2}}) \quad \text{in 2D} \\
    \forall \mathbf{a}\in\Omega, \quad |w-\pi_{\mathcal{M}}w|(\mathbf{a})&=\frac{1}{10}\text{tr}(\mathcal{M}(\mathbf{a})^{-\frac{1}{2}} |H(\mathbf{a})| \mathcal{M}(\mathbf{a})^{-\frac{1}{2}}) \quad \text{in 3D}
\end{aligned}
$$

The continuous linear interpolation model in $L^p$-norm writes:

$$
\mathcal{E}_{L^p}(\mathcal{M}) = \left(\int_{\Omega}|w-\pi_{\mathcal{M}}w|^p d\Omega\right)^{\frac{1}{p}} =  C_d \left(\int_{\Omega}|\text{tr}(\mathcal{M}^{-\frac{1}{2}} |H| \mathcal{M}^{-\frac{1}{2}})|^p d\Omega\right)^{\frac{1}{p}}
$$

Where $C_d$ is a constant depending on the dimension $d$ of the simplices composing $\mathcal{T}_h$.

The solution $\mathcal{M}^{opt}_{L^p}$ of the problem exists and is given by:

$$
\mathcal{M}^{opt}_{L^p} = N^{\frac{2}{d}}\left( \int_{\Omega} \det|H_w|^{\frac{p}{2p+d}} \right)^{-\frac{2}{d}} \det|H_w|^{-\frac{1}{2p+d}}|H_w|
$$

The proofs are given in (Loseille 2011b).

## Hessian computation

In practice, two Hessian-recovery strategies are used in the lib:

- a Least-Squares (LS) approach applied twice (gradient then Hessian), with local systems solved by QR decomposition,
- a Clément-type local $L^2$ projection approach (`L2proj`) (Clement, 1975), similar to the recovery strategy described in Alauzet (2021).

Both approaches operate on vertex values. When the source field is element-centered, a volumetric interpolation to mesh vertices is performed first.