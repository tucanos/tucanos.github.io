# Curvature metric

A metric field can be defined on the boundaries of the computational domain to ensure an accurate representation of the curvature
  - the curvature tensor is computed on an accurate geometrical model, giving principal directions $\mathbf{t_0}$ and $\mathbf{t_1}$ and the associated curvature radii $r_0$ and $r_1$.
  - along direction $\mathbf{t_i}$, a size $h_i = r_i / c$ is imposed (with a given value of $c$, typically $c \sim 5$)
  - along the normal direction $\mathbf n = \mathbf{t_0} \times \mathbf{t_1}$ the edge length $h_n$ can be chosen arbitrarily. It the representation of the curvature is the only objective, one can chose $h_n = \min(h_0, h_1)$.

This metric field, defined on the boundary of the computational domain, can then be extended into the domain using a fixed [gradation](#metric-gradation)

## boundary-normal size for turbulent boundary layers

For the simulation of turbulent boundary layers, it may be interesting to ensure that the mesh will be fine enough to capture the linear (or log if a wall model is available) layer. In this case, the boundary normal size $h_n$ can be chosen to have a target value of $y^+$.

In order to estimate the wall shear stress, required to estimate the target edge length in the normal direction, a wall model should be used (as the initial mesh is in practive much too coarse to resolve the turbulent boundary layer)
