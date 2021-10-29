module Mesh

using DocStringExtensions

export Mesh2D,
    Mesh3D,
    AbstractMesh,
    equispaced_rectangular_mesh,
    sphere_mesh,
    cube_panel_mesh,
    AbstractSphereWarp,
    EquiangularSphereWarp,
    EquidistantSphereWarp,
    unfold_cube_panel_to_plane

abstract type AbstractMesh{I,FT} end

abstract type AbstractSphereWarp end
struct EquiangularSphereWarp <: AbstractSphereWarp end
struct EquidistantSphereWarp <: AbstractSphereWarp end

"""
    Mesh2D{I,IA2D,FT,FTA2D} <: AbstractMesh{FT}

Conformal mesh for a 2D manifold. The manifold can be 
embedded in a higher dimensional space.


                        Quadrilateral

                v4            f3           v3
                  o------------------------o
                  |                        |      face    vertices
                  |       <----------      |             
                  |                 |      |        f1 =>  v1 v2 
               f4 |                 |      | f2     f2 =>  v2 v3
                  |                 |      |        f3 =>  v3 v4
                  |                 |      |        f4 =>  v4 v1
                  |       -----------      |
                  |                        |
                  o------------------------o
                 v1           f1           v2
   
Reference:
https://gsjaardema.github.io/seacas-docs/html/element_types.html#ordering

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Mesh2D{I,IA1D,IA2D,FT,FTA2D,SNT,NB} <: AbstractMesh{I,FT}
    "number of unique vertices in the mesh"
    nverts::I
    "number of unique faces in the mesh"
    nfaces::I
    "number of elements in the mesh"
    nelems::I
    "number of zones in the mesh"
    nbndry::I
    "x₁, x₂, ... coordinates of vertices `(nverts, dim)`, dim can be greater than 2 for 2D manifolds embedded in higher dimensional space"
    coordinates::FTA2D
    "unique vertices `(n_uniquevertices)`"
    unique_verts::IA1D
    "connectivity information for unique vertices"
    uverts_conn::IA1D
    "offset information for uverts_conn `(n_unique_verts + 1)`"
    uverts_offset::IA1D
    "face vertices numbers `(nfaces, 2)`"
    face_verts::IA2D
    "boundary elems for each face `(nfaces, 5)` -> [elem1, localface1, elem2, localface2, relative orientation]"
    face_neighbors::IA2D
    "face numbers on each boundary `(nfaces)`"
    face_boundary::IA1D
    "boundary tags"
    boundary_tags::IA1D
    "boundary tag names"
    boundary_tag_names::SNT
    "face boundary offset for each boundary for face_boundary array"
    face_boundary_offset::IA1D
    "vertices numbers for each elem `(nelems, 4)`"
    elem_verts::IA2D
    "face numbers for each elem `(nelems, 4)`"
    elem_faces::IA2D
end

Mesh2D(
    nverts,
    nfaces,
    nelems,
    nbndry,
    coordinates,
    unique_verts,
    uverts_conn,
    uverts_offset,
    face_verts,
    face_neighbors,
    face_boundary,
    boundary_tags,
    boundary_tag_names,
    face_boundary_offset,
    elem_verts,
    elem_faces,
) = Mesh2D{
    eltype(nverts),
    typeof(face_boundary),
    typeof(face_verts),
    eltype(coordinates),
    typeof(coordinates),
    typeof(boundary_tag_names),
    length(boundary_tag_names),
}(
    nverts,
    nfaces,
    nelems,
    nbndry,
    coordinates,
    unique_verts,
    uverts_conn,
    uverts_offset,
    face_verts,
    face_neighbors,
    face_boundary,
    boundary_tags,
    boundary_tag_names,
    face_boundary_offset,
    elem_verts,
    elem_faces,
)

"""
    Mesh3D{I,IA2D,FT,FTA2D} <: AbstractMesh{FT}

Conformal mesh for 3D manifold.

    
                           Hexahedron
                     v8                     v7      
                      o--------------------o       face no     loop     normal
                     /|                   /|                                   
                    / |                  / |       face 1:=>  1 2 6 5     in (-y)
                   /  |                 /  |       face 2:=>  2 3 7 6    out (+x)
                  /   |                /   |       face 3:=>  3 4 8 7    out (+y)
                 /    |               /    |       face 4:=>  1 5 8 4     in (-x)
                /     |            v6/     |       face 5:=>  1 4 3 2    out (-z)
            v5 o--------------------o      |       face 6:=>  5 6 7 8     in (+z)
               |    v4o-------------|------o v3 
               |     /              |     /        edge no    vertices
               |    /               |    /                              
               |   /                |   /          edge  1:=>  1 2
               |  /                 |  /           edge  2:=>  2 3
               | /                  | /            edge  3:=>  3 4 
               |/                   |/             edge  4:=>  4 1
               o--------------------o              edge  5:=>  1 5
              v1                    v2             edge  6:=>  2 6 
                                                   edge  7:=>  3 7 
                                                   edge  8:=>  4 8
                                                   edge  9:=>  5 6
                                                   edge 10:=>  6 7
                                                   edge 11:=>  7 8
                                                   edge 12:=>  8 5

https://gsjaardema.github.io/seacas-docs/html/element_types.html#ordering

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Mesh3D{I,IA1D,IA2D,FT,FTA2D,SNT,NB} <: AbstractMesh{I,FT}
    "number of nodes in the mesh"
    nverts::I
    "number of edges in the mesh"
    nedges::I
    "number of faces in the mesh"
    nfaces::I
    "number of elements in the mesh"
    nelems::I
    "number of zones in the mesh"
    nbndry::I
    "x₁, x₂, x₃ coordinates of nodes `(nverts, 3)`"
    coordinates::FTA2D
    "unique vertices `(n_uniquevertices)`"
    unique_verts::IA1D
    "connectivity information for unique vertices"
    uverts_conn::IA1D
    "offset information for uverts_conn `(n_unique_verts + 1)`"
    uverts_offset::IA1D
    "edge vertex numbers `(nedges, 2)`"
    edge_verts::IA2D
    "face vertices numbers `(nfaces, 4)`"
    face_verts::IA2D
    "face edge numbers `(nfaces, 4)`"
    face_edges::IA2D
    "boundary elems for each face `(nfaces, 5)` -> [elem1, localface1, elem2, localface2, relative orientation]"
    face_neighbors::IA2D
    "face numbers on each boundary `(nfaces)`"
    face_boundary::IA1D
    "boundary tags"
    boundary_tags::IA1D
    "boundary tag names"
    boundary_tag_names::SNT
    "face boundary offset for each boundary for face_boundary array"
    face_boundary_offset::IA1D
    "vertices numbers for each elem `(nelems, 4)`"
    elem_verts::IA2D
    "face numbers for each elem `(nelems, 4)`"
    elem_faces::IA2D
    "edge numbers for each elem `(nelems, 12)`"
    elem_edges::IA2D
end

Mesh3D(
    nverts,
    nedges,
    nfaces,
    nelems,
    nbndry,
    coordinates,
    unique_verts,
    uverts_conn,
    uverts_offset,
    edge_verts,
    face_verts,
    face_edges,
    face_neighbors,
    face_boundary,
    boundary_tags,
    boundary_tag_names,
    face_boundary_offset,
    elem_verts,
    elem_faces,
    elem_edges,
) = Mesh3D{
    eltype(nverts),
    typeof(face_boundary),
    typeof(face_verts),
    eltype(coordinates),
    typeof(coordinates),
    typeof(boundary_tag_names),
    length(boundary_tag_names),
}(
    nverts,
    nedges,
    nfaces,
    nelems,
    nbndry,
    coordinates,
    unique_verts,
    uverts_conn,
    uverts_offset,
    edge_verts,
    face_verts,
    face_edges,
    face_neighbors,
    face_boundary,
    boundary_tags,
    boundary_tag_names,
    face_boundary_offset,
    elem_verts,
    elem_faces,
    elem_edges,
)

include("box_mesh.jl")
include("sphere_mesh.jl")
include("mesh_utilities.jl")
end
