
# Cavity operators

The remesher is built around local cavities, reconstructed and validated before writing any modification to the global mesh.

## Cavity definition

Given a seed entity $e$ (vertex or edge), the cavity $C(e)$ is the set of elements incident to $e$.

In the implementation:

- `init_from_vertex` builds a cavity from all elements around one vertex,
- `init_from_edge` builds a cavity from elements sharing one edge,
- global vertex IDs are mapped to local IDs (`local2global`) so reconstruction can be handled locally,
- boundary/interface information is split between internal faces, tagged faces, and tagged-face boundaries.

## Cavity reconstruction

For reconstruction, the boundary of the cavity is used to create a new local triangulation/tetrahedralization from a reconstruction vertex $\mathbf v$:

$$
\mathcal F(\mathbf v, C(e)) =
\left\{
K=(\mathbf v,\mathbf g_1,\ldots,\mathbf g_d)
\;\middle|\;
g=(\mathbf g_1,\ldots,\mathbf g_d) \in \partial C(e),\; \mathbf v \notin g
\right\}.
$$


### Filled cavity types

A filled cavity is evaluated with three reconstruction modes:

1. existing cavity vertex (`ExistingVertex`) for collapse and swap,
2. moved cavity vertex (`MovedVertex`) for smoothing,
3. new edge-center vertex (`EdgeCenter`) for split.

The three modes differ in which checks are meaningful: the length, quality and normal checks apply to every reconstruction, while the topology-preservation check only applies when the reconstruction vertex already belongs to the cavity (`ExistingVertex`, i.e. collapse and swap). The checks are described in [Validation of a reconstructed cavity](#validation-of-a-reconstructed-cavity).

## Validation of a reconstructed cavity

A reconstructed cavity is accepted only if its new elements are geometrically valid, close enough to the requested unit size in the metric, and consistent with the tagged geometry. The checks below are predicates evaluated on the filled cavity, before any change is written to the global mesh:

1. metric length and quality bounds (`FilledCavity::check`), see [Metric lengths and quality](#metric-lengths-and-quality),
2. topology preservation (`FilledCavity::check_topo`), see [Topology preservation](#topology-preservation),
3. boundary-normal consistency with the geometry (`FilledCavity::check_normals`), see [Boundary normals](#boundary-normals).

### Metric lengths and quality

For each element $K=(\mathbf v,\mathbf g_1,\ldots,\mathbf g_d)$ of $\mathcal F(\mathbf v,C(e))$, the reconstruction vertex $\mathbf v$ is joined to the vertices $\mathbf x$ of the boundary face $g=(\mathbf g_1,\ldots,\mathbf g_d)$ by new edges. `check(l_min, l_max, q_min)` tests

- the metric edge length of every new edge,
  $$ l = l_{\mathcal M}(\mathbf v,\mathbf x), $$
  against the lower and upper bounds $l_{\min} \le l \le l_{\max}$,
- the metric quality $q(K)$ of every new element against the lower bound $q(K) > q_{\min}$.

The metric edge length uses the metrics at both endpoints, combined by the logarithmic mean of [Numerical computation of lengths](riemannian_metric.md#numerical-computation-of-lengths). The quality is the one of [Element quality](metric.md#element-quality), computed with the minimum metric over the vertices of $K$. Only the edges and elements created by the reconstruction are tested: the boundary faces of the cavity are faces of existing elements and are not re-evaluated.

The returned status is determined by the first failing edge or element:

| Condition | Status | Meaning |
| --- | --- | --- |
| $l < l_{\min}$ | `ShortEdge(l)` | the reconstruction would create an edge shorter than requested |
| $l > l_{\max}$ | `LongEdge(l)` | the reconstruction would create an edge longer than requested |
| $q(K) < 0$ | `Invalid` | the element has a non-positive measure (inverted or flat), see [Element quality](metric.md#element-quality) |
| $q(K) \le q_{\min}$ | `LowQuality(q)` | the element is valid but not good enough |
| otherwise | `Ok(q_{\min}^{\text{new}})` | all elements are valid; $q_{\min}^{\text{new}}$ is their minimum quality |

The bounds are not symmetric in their strictness: $l = l_{\min}$ and $l = l_{\max}$ pass, while $q = q_{\min}$ fails. The minimum quality reported by `Ok` is used to rank swap candidates and to update the threshold inside a swap loop, see [Per-operator thresholds](#per-operator-thresholds).

### Topology preservation

The reconstruction replaces the seed $e$ by $\mathbf v$ and re-creates the elements around it. `check_topo` ensures that this local rewrite does not duplicate an entity nor make an edge incident to more entities than before. It is only defined for `ExistingVertex`, because smoothing keeps the connectivity unchanged and splitting inserts a vertex that belongs to no existing element: neither can duplicate a topological entity. For collapse and swap, the checked entity depends on the dimension of the elements:

| Mesh | Duplicate entity | Pinched entity |
| --- | --- | --- |
| volume ($d = D$, tetrahedra) | tagged face | edge shared by tagged faces (collapse only) |
| surface ($d = D-1$, triangles) | element | edge shared by elements |

- **Duplicate entity.** In a surface mesh, the element $(\mathbf v,g)$ built from $\mathbf v$ and a boundary face $g$ is rejected if it is already incident to $\mathbf v$, which would create two identical triangles. In a volume mesh, the face $(\mathbf v,b)$ built from $\mathbf v$ and an edge $b$ of a tagged cavity boundary face is rejected if that sorted face already belongs to the tagged-face map, which would create two coincident tagged patches.
- **Pinched edge, collapse.** Collapsing $j$ onto $i$ merges the edge $(j,k)$ into $(i,k)$ for every boundary vertex $k$. Counting, before the collapse, the tagged faces (volume) or elements (surface) that contain $(j,k)$, and after the collapse those that contain $(i,k)$ or $(j,k)$ but not both $i$ and $j$ (the latter belong to removed elements), the reconstruction is rejected when
  $$ \mathrm{card}_{\text{after}} > \mathrm{card}_{\text{before}}. $$
  This prevents an edge from becoming incident to more tagged faces, or to more elements, than before the collapse.
- **Pinched edge, swap.** Swapping the edge of a triangle cavity replaces it by an edge $(\mathbf v,k)$ for each cavity boundary vertex $k$. If $(\mathbf v,k)$ is already an edge of the mesh, the number of elements incident to it after the swap is the number of reconstructed elements that contain it plus the number of existing elements outside the cavity that contain it; the candidate is rejected when this sum exceeds two, since a surface edge may only be shared by two triangles. If $(\mathbf v,k)$ is not an existing edge, no pinch can occur and the test passes. In a volume mesh this test is not applied to swaps, because an interior edge is legitimately shared by any number of tetrahedra; only the duplicate tagged-face test remains.

`check_topo` only handles duplication and pinching. Tag compatibility (which vertex may replace which entity, and which tags must be preserved) is checked by the operators before the cavity is built, using the topology hierarchy of [Locating an edge from tags](topology.md#locating-an-edge-from-tags).

> [!TIP]
> **Cavity construction and topology check messages**
>
> | Message | Meaning / action |
> | --- | --- |
> | `Empty vertex cavity: i = {i:?} (n_verts = {} / n_elems = {})` | `init_from_vertex` was called on a vertex referenced by no element. Do not seed a cavity from an isolated vertex. |
> | `Empty edge cavity: edg = {edg:?} ({})--> elements = {:?} & {:?}` | `init_from_edge` was called on an edge shared by no element. The edge is not part of the remeshed domain. |
> | `not implemented` | `check_topo` was called on elements whose dimension is neither $d$ (volume) nor $d-1$ (surface); only these two cases are defined. |
> | `internal error: entered unreachable code` | `check_topo` was called with a `MovedVertex` or `EdgeCenter` filled cavity. Only collapse and swap have a reconstruction vertex that already belongs to the cavity. |
> | `assertion failed: ...` in a pinch helper | Internal invariant: the edge joining the reconstruction vertex to a cavity boundary vertex, or its incident elements or tagged faces, is missing from the remesher. This indicates a corrupt state; report a bug. |

### Boundary normals

A reconstruction that is valid in the metric can still create a boundary that deviates from the tagged geometry. `check_normals(topo, geom, max_angle)` compares the normal of each new boundary face with the normal of the tagged patch it belongs to, and rejects the reconstruction if the angle between them exceeds `max_angle` (in degrees):

- **surface mesh** ($d = D-1$): every reconstructed element is a new face of the surface; its center and normal are compared with the geometry normal at the tag of its boundary face,
- **volume mesh** ($d = D$): only the faces of the reconstructed elements that lie on a tagged boundary are checked, i.e. the faces spanned by $\mathbf v$ and an edge of a tagged cavity boundary face. The orientation of the cavity boundary face is used to flip the normal when needed, so the comparison is always made with the outward normal of the new element on that patch.

The angle is computed as $\arccos(\mathbf n\cdot\mathbf n_{\text{ref}})$ in degrees, with both normals normalized, and is available from the geometry of [Meshed geometries](geometry.md#meshed-geometries). This is the same `max_angle` parameter for all four operators (default 25 degrees); it is checked after the topology and before the length/quality ranking in swap.

> [!TIP]
> **Boundary-normal check message**
>
> | Message | Meaning / action |
> | --- | --- |
> | `Invalid tag {}` | The volume branch of `check_normals` found a tagged cavity face with a non-positive tag. Negative tags mark fixed entities and must be rejected by the operator before the cavity is built; a zero tag is not a valid face tag. |
>
> The call to the geometry can also raise the [`Invalid face tag`](geometry.md#meshed-geometries) error when no patch is registered for the tag of a new boundary face.

### Per-operator thresholds

The bounds passed to `check` are derived from the current cavity, so that a reconstruction is compared with the elements it replaces. Let

- $L_{\min}$ and $L_{\max}$ be the minimum and maximum metric edge lengths over the edges of the current cavity elements,
- $Q_{\min}$ be the minimum metric quality over the current cavity elements,

then the thresholds are:

- **Collapse**: $l_{\min}=0$, $l_{\max}=\max(\texttt{max\_l\_abs},\,\texttt{max\_l\_rel}\,L_{\max})$, $q_{\min}=\min(\texttt{min\_q\_abs},\,\texttt{min\_q\_rel}\,Q_{\min})$.
- **Split**: $l_{\min}=\min(\texttt{min\_l\_abs},\,\texttt{min\_l\_rel}\,L_{\min})$, $l_{\max}=+\infty$, $q_{\min}=\min(\texttt{min\_q\_abs},\,r\,Q_{\min})$.
- **Swap**: $l_{\min}=\min(\texttt{min\_l\_abs},\,\texttt{min\_l\_rel}\,L_{\min})$, $l_{\max}=\max(\texttt{max\_l\_abs},\,\texttt{max\_l\_rel}\,L_{\max})$, $q_{\min}=Q_{\min}$, raised to the minimum quality of each accepted candidate.
- **Smooth**: $l_{\min}=0$, $l_{\max}=+\infty$, $q_{\min}=Q_{\min}$.

For split, the relative quality factor $r$ is `min_q_rel` when the split edge is interior (`edge_tag` of dimension $d$) and `min_q_rel_bdy` when it lies on a tagged boundary. Collapse and split therefore tolerate a quality or length degradation bounded by the relative factors, whereas swap and smooth require the new minimum quality to be strictly above $Q_{\min}$; for swap, raising $q_{\min}$ to the quality of the last accepted candidate makes the loop select the candidate with the largest new minimum quality.

Default parameter sets use `max_l_rel = 1` for collapse and `1.5` for swap, `max_l_abs` $= 1.5\sqrt{2}$, `min_l_rel = 1` for split and `0.75` for swap, `min_l_abs` $= 0.75/\sqrt{2}$, `min_q_rel = 1` for collapse and `0.8` for split, `min_q_rel_bdy = 0.5`, `min_q_abs = 0.5` for collapse and `0.3` for split, and `max_angle = 25` degrees for all operators.

### Failure handling

A rejected candidate leaves the cavity and the global mesh unchanged, and the operator counts the failure. `ShortEdge`, `LongEdge` and `LowQuality` reject the split directly. Only an `Invalid` split on a boundary edge (`edge_tag` of dimension $< d$) is retried on a larger cavity, where the [cavity extension](#cavity-extension) and [element suppression](#element-suppression) mechanisms attempt to remove the inverted configuration. In debug mode, the whole remesher state is re-validated by `Remesher::check` at the end of a step; its messages are listed in [Validity](validity.md#remesher-consistency-errors).

## Operator-specific cavity usage

| Operation | Seed | Reconstruction | Example (2D) |
| --- | --- | --- | --- |
| Collapse | removed endpoint| kept endpoint | ![collapse](images/collapse-cavity.svg) |
| Split| edge | (projected) edge center | ![split](images/split-cavity.svg) |
| Swap| edge | the cavity boundary vertex that maximizes post-swap minimum quality. | ![swap](images/swap-cavity.svg) |


## Cavity extension 

When splitting boundary edges, it may happen that reconstructing the cavity from the projection of the edge center would create invalid elements. In this case `Cavity::extend` grows the cavity across a selected cavity boundary face by adding the adjacent element on the other side (same element tag required). 

![extension](images/extension-cavity.svg) 

Extension is a split operation and is repeated up to `SplitParams::max_extensions` times (default 20) until all candidate reconstructed elements are valid.

> [!TIP]
> **Cavity extension message**
>
> | Message | Meaning / action |
> | --- | --- |
> | `face {f:?} belongs to {} elements` | `extend` tried to cross a cavity boundary face shared by more than two elements. The remesher assumes manifold faces; a non-manifold face cannot be crossed by extension (see [Non-manifold meshes](non_manifold.md)). |

## Element suppression

> [!WARNING]
> Element suppression is currently implemented for **3D (tetrahedral) meshes only**. It is disabled in 2D.

A tetrahedron is removed when its cavity contains exactly two tagged faces that share the same positive tag, the other two cavity faces are untagged, and re-tagging those two opposite faces with that tag improves the geometric representation. Geometrically, the largest normal deviation `max(a_0, a_1)` of the two tagged faces must exceed the largest deviation `max(a_2, a_3)` of the two untagged faces (angles are measured against the geometry normals for the tag). The two original tagged faces are then removed, the two opposite faces are tagged, and the element is deleted.

![extension](images/extension-cavity-suppression.svg) 

