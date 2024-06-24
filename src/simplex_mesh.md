# Simplex meshes

Mesh adaptation can only be applied to conformal meshes made only of simplex elements, i.e. triangles in 2D and tetrahedra in 3d. For other types of meshes, the elements must initially be split into simplices: algorithms for standard elements are given in [*How to Subdivide Pyramids, Prisms and Hexahedra into Tetrahedra*, Julien Dompierre Paul Labbé Marie-Gabrielle Vallet Ricardo Camarero](https://www.researchgate.net/publication/221561839_How_to_Subdivide_Pyramids_Prisms_and_Hexahedra_into_Tetrahedra) and implemented in `tucanos`.


## Mesh

A mesh consists of
- Vertices $\mathbf x_i \in \mathbb R^d$
- Elements $K_i$ and element tags $T_{K, i}$
- Tagged faces $F_i$ and face tags $T_{F, i}$

It is characterized by a spatial dimension $d$ as well as a cell dimension $d_{cell} \le d$ (3 for tetrahedra, 2 for triangles, ...).
The set of mesh entities of dimension $d$ is denoted $\mathcal E(d)$.

A mesh is valid if
- 
- All the elements have a volume $|K_i| > 0$,
- All the faces must be either either
    - tagged if 
        - belong to only one element. In this case, they must be oriented outwards
        - belong to two elements with different tags
        - (belong to more than two elements)
    - untagged

and consistent with the geometry (i.e. a CAD-like model) if all the vertices of the tagged faces lie on the corresponding geometry surfaces.
