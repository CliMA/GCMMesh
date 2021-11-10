# utilities for Mesh2D
function opposing_face(mesh::Mesh2D, elem, local_face)
    @assert 0 < elem ≤ mesh.nelems
    @assert 0 < local_face ≤ 4

    neighbors = mesh.face_neighbors
    face = mesh.elem_faces[elem, local_face]

    if neighbors[face, 1] == elem && neighbors[face, 2] == local_face
        neigh_elem, neigh_lface = neighbors[face, 3], neighbors[face, 4]
    elseif neighbors[face, 3] == elem && neighbors[face, 4] == local_face
        neigh_elem, neigh_lface = neighbors[face, 1], neighbors[face, 2]
    else
        error(
            "opposing_face: Fatal error, elem not found in face_neighbors; face_neighbors[$face] = $(neighbors[face, :])",
        )
    end
    if neigh_elem == 0
        neigh_lface = local_face
    end
    reversed = neighbors[face, 5] == -1
    return neigh_elem, neigh_lface, reversed
end

function interior_faces(mesh::Mesh2D)
    tags = mesh.boundary_tags
    offset = mesh.face_boundary_offset

    if tags[1] ≠ 0
        return nothing
    else
        return mesh.face_boundary[1:offset[2]-1]
    end
end

function boundary_faces(mesh::Mesh2D, btag::Int)
    tags = mesh.boundary_tags
    offset = mesh.face_boundary_offset
    loc = findfirst(tags .== btag)
    if loc == nothing
        return nothing
    else
        return mesh.face_boundary[offset[loc]:offset[loc+1]-1]
    end
end

function boundary_face_info(mesh::Mesh2D, face::Int)
    face_neighbors = mesh.face_neighbors
    if face_neighbors[face, 1] == 0
        return tuple(face_neighbors[face, 3:4]...)
    else
        return tuple(face_neighbors[face, 1:2]...)
    end
end

function boundary_tag(mesh::Mesh2D, tag_name::Symbol)
    tag_names = mesh.boundary_tag_names
    tags = mesh.boundary_tags
    loc = findfirst(mesh.boundary_tag_names .== tag_name)
    if loc == nothing
        return nothing
    else
        return tags[loc]
    end
end

function unique_vertex_connectivity(mesh::Mesh2D, uvertno::Int)
    offset = mesh.uverts_offset
    conn = mesh.uverts_conn
    st, en = offset[uvertno], offset[uvertno+1]
    nconn = div(en - st, 2)

    if nconn == 0
        return nothing
    else
        return [tuple(conn[(st+(i-1)*2):(st+(i-1)*2+1)]...) for i = 1:nconn]
    end
end

function vertex_coordinates(mesh::Mesh2D, elem)
    @assert 0 < elem ≤ mesh.nelems
    verts = mesh.elem_verts[elem, :]
    coords = mesh.coordinates
    return [tuple(coords[vert, :]...) for vert in verts]
end

# utilities for Mesh1D
# Note: for 1D meshes, faces are identical to vertices
function opposing_face(mesh::Mesh1D, elem, local_face)
    @assert 0 < elem ≤ mesh.nelems
    @assert 0 < local_face ≤ 2

    neighbors = mesh.vertex_neighbors
    face = mesh.elem_verts[elem, local_face]

    if neighbors[face, 1] == elem && neighbors[face, 2] == local_face
        neigh_elem, neigh_lface = neighbors[face, 3], neighbors[face, 4]
    elseif neighbors[face, 3] == elem && neighbors[face, 4] == local_face
        neigh_elem, neigh_lface = neighbors[face, 1], neighbors[face, 2]
    else
        error(
            "opposing_face: Fatal error, elem not found in face_neighbors; face_neighbors[$face] = $(neighbors[face, :])",
        )
    end
    if neigh_elem == 0
        neigh_lface = local_face
    end
    return neigh_elem, neigh_lface
end

function interior_faces(mesh::Mesh1D)
    tags = mesh.boundary_tags
    offset = mesh.vertex_boundary_offset

    if tags[1] ≠ 0
        return nothing
    else
        return mesh.vertex_boundary[1:offset[2]-1]
    end
end

function boundary_faces(mesh::Mesh1D, btag::Int)
    tags = mesh.boundary_tags
    offset = mesh.vertex_boundary_offset
    loc = findfirst(tags .== btag)
    if loc == nothing
        return nothing
    else
        return mesh.vertex_boundary[offset[loc]:offset[loc+1]-1]
    end
end

function boundary_face_info(mesh::Mesh1D, face::Int)
    face_neighbors = mesh.vertex_neighbors
    if face_neighbors[face, 1] == 0
        return tuple(face_neighbors[face, 3:4]...)
    else
        return tuple(face_neighbors[face, 1:2]...)
    end
end

function boundary_tag(mesh::Mesh1D, tag_name::Symbol)
    tag_names = mesh.boundary_tag_names
    tags = mesh.boundary_tags
    loc = findfirst(mesh.boundary_tag_names .== tag_name)
    if loc == nothing
        return nothing
    else
        return tags[loc]
    end
end

function unique_vertex_connectivity(mesh::Mesh1D, uvertno::Int)
    neighbors = mesh.vertex_neighbors
    vt = mesh.unique_verts[uvertno]
    if neighbors[vt, 1] == 0 && neighbors[vt, 3] ≠ 0
        return [tuple(neighbors[vt, 3:4]...)]
    elseif neighbors[vt, 1] ≠ 0 && neighbors[vt, 3] == 0
        return [tuple(neighbors[vt, 1:2]...)]
    else
        return [tuple(neighbors[vt, 1:2]...), tuple(neighbors[vt, 3:4]...)]
    end
end

function vertex_coordinates(mesh::Mesh1D, elem)
    @assert 0 < elem ≤ mesh.nelems
    verts = mesh.elem_verts[elem, :]
    coords = mesh.coordinates
    return [tuple(coords[vert, :]...) for vert in verts]
end
