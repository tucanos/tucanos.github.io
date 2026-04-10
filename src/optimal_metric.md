# Optimal metric definition
## Feature based adaptation
Feature-based in metric-based adaptation consists in applying the **multiscale error estimator** in $L^p$-norm to a scalar function in the domain. This function is denoted as a sensor, or a feature, on which the mesh is adapted to reduce the global interpolation error on this particular sensor, while automatically capturing all scales. In CFD computations, the Mach number is typically used because it captures features of interest such as boundary layers, shocks, and wakes.

Let $\mathcal{T}_h$ be a simplex mesh composed of simplex elements $K$. Let $E_{L^p}(\mathcal{T}_h)$ be an error model in $L^p$-norm. The optimal mesh $\mathcal{T}_h^{opt}$ is given by the following minimization problem with constraint $\mathcal{C}(\mathcal{T}_h)=N$:

$$
\mathcal{T}_h^{opt}=\operatorname*{argmin}_{\mathcal{C}(\mathcal{T}_h)=N} {E_{L^p}(\mathcal{T}_h)}
$$

In the continuous mesh theoretical framework, the goal is to find the optimal continuous mesh $\mathcal{M}_{opt}$ which minimizes the given continuous error model $\mathcal{E}_{L^p}$ for a fixed continuous mesh complexity $C(\mathcal{M}) = N$:

$$
\mathcal{M}_{opt}=\operatorname*{argmin}_{\mathcal{C}(\mathcal{M})=N} {\mathcal{E}_{L^p}(\mathcal{M})}
$$

## Multiscale error estimator in $L^p$-norm and feature-based mesh adaptation

Let $w=f(\mathbf{u})$ be a twice derivable smooth function of the state $\mathbf{u}$. Let $\mathcal{T}_h$ be a simplex mesh on which the function $w$ is piecewise linear approximated on elements $K$. The piecewise approximation is defined by the linear interpolate $\Pi_h w$ on the mesh $\mathcal{T}_h$. Then, the global interpolation error of $w$ in $\mathcal{T}_h$ is defined by the $L^p$-norm of its local interpolation errors, namely:

$$
E_{L^p}(\mathcal{T}_h) = \left(\int_{\Omega}|w-\Pi_hw|^p d\Omega\right)^{\frac{1}{p}}.
$$

In (Loseille 2011), the authors prove the existence and unicity of the continuous linear interpolate $\pi_{\mathcal{M}}w$ that defines the continuous linear interpolation error related to the continuous mesh $(\mathcal{M}(\mathbf{x}))_{\mathbf{x}\in\Omega}$. From corollary 3.4 of (Loseille 2011), the following continuous linear interpolation estimate holds in 2D and 3D:

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
Tucanos features two methods to compute the Hessian of a given scalar quantity. The first one uses a simple Least-Squares (LS) formulation that is performed twice. The resolution of the Least-squares is carried using a QR decomposition. The second method uses Clément’s L2-projection local operator (L2proj) (Clement, 1975) similarly to the Hessian recovery method presented in (Alauzet, 2021). These methods use the values at the vertices of the mesh, obtained via a volumetric interpolation of the values from the elements to the nodes. This interpolation is performed using the Tucanos library.