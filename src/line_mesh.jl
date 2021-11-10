"""
    equispaced_line_mesh(
    x1min::FT,
    x1max::FT,
    n1,
    per = false,
) where {FT<:AbstractFloat}

This function generates an equispaced line mesh with user provided coordinates for left and right
domain vertices, number of mesh elements and periodicity information.
"""
function equispaced_line_mesh(
    x1min::FT,
    x1max::FT,
    n1,
    per = false,
) where {FT<:AbstractFloat}
    x1c = range(x1min, x1max; length = n1 + 1)
    return Mesh1D(x1c, per)
end

function Mesh1D(x1c, per = false)
    I = Int
    FT = eltype(x1c)
    nverts = length(x1c)
    nelems = nverts - 1
    nbndry = per ? 1 : 3
    coordinates = Array{FT}(undef, nverts, 1)
    coordinates[:, 1] .= x1c
    vertex_neighbors = zeros(I, nverts, 4)
    tag_names = (:interior, :left, :right)
    if per
        unique_verts = Array{I}(1:nverts-1)
        boundary_tag_names = tuple(tag_names[1])
        boundary_tags = [0]
        vertex_neighbors[1, :] .= [nelems, 2, 1, 1]
        vertex_neighbors[nverts, :] .= [nelems, 2, 1, 1]
        vertex_boundary = Vector(1:nverts-1)
        vertex_boundary_offset = [1, nverts]
    else
        unique_verts = Array{I}(1:nverts)
        vertex_neighbors[1, :] .= [0, 0, 1, 1]
        vertex_neighbors[nverts, :] .= [nelems, 2, 0, 0]
        if nelems > 1
            boundary_tag_names = tuple(tag_names...)
            boundary_tags = [0, 1, 2]
            vertex_boundary = vcat(Vector(2:nverts-1), 1, nverts)
            vertex_boundary_offset = [1, nverts - 1, nverts, nverts + 1]
        else
            boundary_tag_names = tuple(tag_names[2:end]...)
            boundary_tags = [1, 2]
            vertex_boundary = vcat(1, nverts)
            vertex_boundary_offset = [1, nverts, nverts + 1]
        end
    end
    for vt = 2:nverts-1
        vertex_neighbors[vt, :] .= [vt - 1, 2, vt, 1]
    end
    elem_verts = Array{I}(undef, nelems, 2)
    for el = 1:nelems
        elem_verts[el, :] .= [el, el + 1]
    end
    return Mesh1D(
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
end
