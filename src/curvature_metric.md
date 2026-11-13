# Curvature metric

A metric field can be defined on boundary patches to control geometric error induced by curvature.

- The curvature tensor is computed on a geometric model, giving principal directions $\mathbf{t_0}$ and $\mathbf{t_1}$ and associated curvature radii $r_0$ and $r_1$.
- Along principal direction $\mathbf{t_i}$, target size is imposed as $h_i = r_i / r_h$, where $r_h$ is the curvature-radius-to-element-size ratio (typically $r_h \approx 4$).
- Along the surface-normal direction $\mathbf n = \mathbf{t_0} \times \mathbf{t_1}$, the size $h_n$ can be chosen independently. If only curvature representation is targeted, a common choice is $h_n = \min(h_0, h_1)$ (this is the default when no normal size is supplied).

This boundary metric is then extended into the volume with fixed [metric gradation](metric_operations.md#metric-gradation), as implemented in `tucanos/src/metric/field/curvature.rs`.

## Boundary-normal size for turbulent boundary layers

For turbulent boundary-layer simulations, it is often necessary to constrain the boundary-normal size to meet a target $y^+$.

In this setting, $h_n$ is chosen from a wall-shear estimate. A wall model is typically required because the initial mesh is usually too coarse to resolve the near-wall layer directly.
