module Mesh

using DocStringExtensions

export Mesh1D,
    Mesh2D,
    Mesh3D,
    AbstractMesh,
    equispaced_line_mesh,
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
    Mesh1D{I,IA1D,IA2D,FT,FTA2D,SNT,NB} <: AbstractMesh{I,FT}

Conformal mesh for a 1D manifold. The manifold can be embedded in a
higher dimensional space.

For 1D meshes, faces are identical to meshes.
"""
struct Mesh1D{I,IA1D,IA2D,FT,FTA2D,SNT,NB} <: AbstractMesh{I,FT}
    "number of unique vertices in the mesh"
    nverts::I
    "number of elements in the mesh"
    nelems::I
    "number of zones in the mesh"
    nbndry::I
    "x₁, x₂, ... coordinates of vertices `(nverts, dim)`, dim can be greater than 1 for 1D manifolds embedded in higher dimensional space"
    coordinates::FTA2D
    "unique vertices `(n_uniquevertices)`"
    unique_verts::IA1D
    "vertex neighbors"
    vertex_neighbors::IA2D
    "vertex numbers on each boundary"
    vertex_boundary::IA1D
    "boundary tags"
    boundary_tags::IA1D
    "boundary tag names"
    boundary_tag_names::SNT
    "vertex boundary offset"
    vertex_boundary_offset::IA1D
    "vertices numbers for each elem `(nelems, 2)`"
    elem_verts::IA2D
end

Mesh1D(
    nverts,
    nelems,
    nbndry,
    coordinates,
    unique_verts,
    vertex_neighbors,
    vertex_boundary,
    boundary_tags,
    boundary_tag_names,
    vertex_boundary_offset,
    elem_verts,
) = Mesh1D{
    eltype(nverts),
    typeof(unique_verts),
    typeof(elem_verts),
    eltype(coordinates),
    typeof(coordinates),
    typeof(boundary_tag_names),
    length(boundary_tag_names),
}(
    nverts,
    nelems,
    nbndry,
    coordinates,
    unique_verts,
    vertex_neighbors,
    vertex_boundary,
    boundary_tags,
    boundary_tag_names,
    vertex_boundary_offset,
    elem_verts,
)
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

include("line_mesh.jl")
include("box_mesh.jl")
include("sphere_mesh.jl")
include("mesh_utilities.jl")
end
