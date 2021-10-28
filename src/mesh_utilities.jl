function opposing_face(mesh::AbstractMesh, elem, local_face)
    @assert 0 < elem ≤ mesh.nelems
    if mesh isa Mesh2D
        @assert 0 < local_face ≤ 4
    elseif mesh isa Mesh3D
        @assert 0 < local_face ≤ 6
    else
        error("opposing face: unknown mesh type $(typeof(mesh))")
    end

    neighbors = mesh.face_neighbors
    face = mesh.elem_faces[elem, local_face]

    if neighbors[face, 1] == elem
        neigh_elem, neigh_lface = neighbors[face, 3], neighbors[face, 4]
    elseif neighbors[face, 3] == elem
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

function boundary_faces(mesh::AbstractMesh, btag::Int)
    tags = mesh.boundary_tags
    offset = mesh.face_boundary_offset
    loc = findfirst(tags .== btag)
    if loc == nothing
        return nothing
    else
        return mesh.face_boundary[offset[loc]:offset[loc+1]-1]
    end
end

function boundary_face_info(mesh::AbstractMesh, face::Int)
    face_neighbors = mesh.face_neighbors
    if face_neighbors[face, 1] == 0
        return tuple(face_neighbors[face, 3:4]...)
    else
        return tuple(face_neighbors[face, 1:2]...)
    end
end

function boundary_tag(mesh::AbstractMesh, tag_name::Symbol)
    tag_names = mesh.boundary_tag_names
    tags = mesh.boundary_tags
    loc = findfirst(mesh.boundary_tag_names .== tag_name)
    if loc == nothing
        return nothing
    else
        return tags[loc]
    end
end

function unique_vertex_connectivity(mesh, uvertno::Int)
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

function vertex_coordinates(mesh::AbstractMesh, elem)
    @assert 0 < elem ≤ mesh.nelems
    verts = mesh.elem_verts[elem, :]
    coords = mesh.coordinates
    return [tuple(coords[vert, :]...) for vert in verts]
end
