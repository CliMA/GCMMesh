"""
    equispaced_rectangular_mesh(
    x1min::FT,
    x2min::FT,
    x1max::FT,
    x2max::FT,
    n1,
    n2,
    per = (false, false),
) where {FT<:AbstractFloat}

This function generates an equispaced rectangular mesh with user provided coordinates of left bottom and
right top vertices, number of mesh elements in each direction and periodicity information.
"""
function equispaced_rectangular_mesh(
    x1min::FT,
    x2min::FT,
    x1max::FT,
    x2max::FT,
    n1,
    n2,
    per = (false, false),
) where {FT<:AbstractFloat}
    x1c = range(x1min, x1max; length = n1 + 1)
    x2c = range(x2min, x2max; length = n2 + 1)
    return Mesh2D(x1c, x2c, per)
end

"""
    Mesh2D(x1c, x2c, per = (false, false))

This function builds a 2D rectangular mesh with points located at x1c and x2c in x1 and x2 directions respectively.
"""
function Mesh2D(x1c, x2c, per = (false, false))
    FT = eltype(x1c)
    I = Int

    nx1, nx2 = length(x1c), length(x2c)
    nel1, nel2 = nx1 - 1, nx2 - 1

    nverts = nx1 * nx2
    nfaces = nx1 * (nx2 - 1) + (nx1 - 1) * nx2
    nelems = (nx1 - 1) * (nx2 - 1)
    nbndry = 4

    nfc1 = nx1 * (nx2 - 1) # faces with normals along x1 direction 

    ndmat = reshape(1:nverts, nx1, nx2) # node numbers

    fcmat1 = reshape(1:nfc1, nx1, nx2 - 1)        # faces with normals along x1 direction
    fcmat2 = reshape((nfc1+1):nfaces, nx1 - 1, nx2) # faces with normals along x2 direction

    emat = reshape(1:nelems, nx1 - 1, nx2 - 1) # element numbers

    coordinates = [repeat(x1c, 1, nx2)[:] repeat(x2c', nx1, 1)[:]] # node coordinates

    face_verts = hcat(
        vcat(ndmat[:, 1:(nx2-1)][:], ndmat[1:(nx1-1), :][:]), # face nodes
        vcat(ndmat[:, 2:nx2][:], ndmat[2:nx1, :][:]),
    )

    face_boundary = Vector{I}(zeros(nfaces))
    # accounting for periodic boundaries, wherever applicable
    periodic = per[1] || per[2]
    periodic_faces = Dict{Int,Int}()
    global_vert = Vector{I}(1:nverts)
    if per[1]
        bdy1 = emat[end:end, :]
        bdy2 = emat[1:1, :]
        nbndry -= 2
        for pfc = 1:(nx2-1)
            pfc1, pfc2 = fcmat1[1, pfc], fcmat1[nx1, pfc]
            periodic_faces[pfc1] = fcmat1[pfc2]
            periodic_faces[pfc2] = fcmat1[pfc1]
        end
        for i = 1:nx2
            global_vert[ndmat[nx1, i]] = ndmat[1, i]
        end
    else
        bdy1 = bdy2 = zeros(I, 1, nx2 - 1)
        face_boundary[fcmat1[1, :]] .= 1  # left boundary
        face_boundary[fcmat1[end, :]] .= 2  # right boundary
    end
    if per[2]
        bdy3 = emat[:, end:end]
        bdy4 = emat[:, 1:1]
        nbndry -= 2
        for pfc = 1:(nx1-1)
            periodic_faces[fcmat2[pfc, 1]] = fcmat2[pfc, nx2]
            periodic_faces[fcmat2[pfc, nx2]] = fcmat2[pfc, 1]
        end
        for i = 1:nx1
            global_vert[ndmat[i, nx2]] = global_vert[ndmat[i, 1]]
        end
    else
        bdy3 = bdy4 = zeros(I, nx1 - 1, 1)
        face_boundary[fcmat2[:, 1]] .= 3  # bottom boundary
        face_boundary[fcmat2[:, end]] .= 4  # top boundary
    end

    face_neighbors = hcat(
        vcat(vcat(bdy1, emat)[:], hcat(bdy3, emat)[:]),
        zeros(I, nfaces),
        vcat(vcat(emat, bdy2)[:], hcat(emat, bdy4)[:]),
        zeros(I, nfaces),
        -ones(I, nfaces), # relative orientation is always reversed
    )

    elem_verts = hcat(
        ndmat[1:(nx1-1), 1:(nx2-1)][:], # node numbers (local vertex # 1)
        ndmat[2:nx1, 1:(nx2-1)][:], # for each element (local vertex # 2)
        ndmat[2:nx1, 2:nx2][:], #(local vertex # 3)
        ndmat[1:(nx1-1), 2:nx2][:], # (local vertex # 4)
    )
    elem_faces = hcat( # face numbers for each element
        fcmat2[:, 1:(nx2-1)][:], # local face 1
        fcmat1[2:nx1, :][:], # local face 2
        fcmat2[:, 2:nx2][:], # local face face 3
        fcmat1[1:(nx1-1), :][:], # local face 4
    )
    ref_fc_verts = [
        1 2 3 4
        2 3 4 1
    ]
    for fc = 1:nfaces
        elems = (face_neighbors[fc, 1], face_neighbors[fc, 3])

        if periodic && haskey(periodic_faces, fc)
            shadow = periodic_faces[fc]
            localface1 = findfirst(elem_faces[elems[1], :] .== fc)
            if !isnothing(localface1)
                localface2 = findfirst(elem_faces[elems[2], :] .== shadow)
                if !isnothing(localface2)
                    face_neighbors[fc, 2] = localface1
                    face_neighbors[fc, 4] = localface2
                else
                    error(
                        "rectangular_mesh: Fatal error, face could not be located in neighboring element",
                    )
                end
            else
                localface1 = findfirst(elem_faces[elems[1], :] .== shadow)
                localface2 = findfirst(elem_faces[elems[2], :] .== fc)
                if !isnothing(localface1) && !isnothing(localface2)
                    face_neighbors[fc, 2] = localface1
                    face_neighbors[fc, 4] = localface2
                else
                    error(
                        "rectangular_mesh: Fatal error, face could not be located in neighboring element",
                    )
                end
            end
        else # no periodic boundaries
            for e = 1:2
                el = elems[e]
                if el ≠ 0
                    localface = findfirst(elem_faces[el, :] .== fc)
                    face_neighbors[fc, 2+(e-1)*2] = localface
                end
            end
        end
    end

    if periodic # marking shadow faces with -1
        for (per, shd) in periodic_faces
            face_boundary[max(per, shd)] = -1
        end
    end
    boundary_tags = sort(unique(face_boundary))
    face_boundary_offset = ones(I, length(boundary_tags) + 1)
    face_boundary_tags = sort(face_boundary)
    face_boundary = sortperm(face_boundary)

    if periodic # remove redundant shadow faces
        loc = findlast(face_boundary_tags .== -1) + 1
        boundary_tags = boundary_tags[2:end]
        face_boundary_offset = face_boundary_offset[2:end]
        face_boundary_tags = face_boundary_tags[loc:end]
        face_boundary = face_boundary[loc:end]
    end

    for i = 1:(length(boundary_tags)-1)
        tag = boundary_tags[i+1]
        face_boundary_offset[i+1] = findfirst(face_boundary_tags .== tag)
    end
    face_boundary_offset[end] = length(face_boundary) + 1
    tag_names = (:interior, :west, :east, :south, :north)

    boundary_tag_names = tuple(tag_names[boundary_tags.+1]...)
    NB = length(boundary_tag_names)

    # add unique vertex iterator information
    vtconn = map(i -> zeros(I, i), zeros(I, nverts))

    for el = 1:nelems
        for lv = 1:4
            vt = global_vert[elem_verts[el, lv]]
            push!(vtconn[vt], el, lv)
        end
    end
    unique_verts = I[]
    uverts_conn = I[]
    uverts_offset = I.([1])
    for vt = 1:nverts
        lconn = length(vtconn[vt])
        if lconn > 0
            push!(unique_verts, vt)
            push!(uverts_conn, vtconn[vt]...)
            push!(uverts_offset, uverts_offset[end] + lconn)
        end
    end

    return Mesh2D(
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
end
