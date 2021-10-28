using Plots

using GCMMesh
using GCMMesh.Mesh
using GCMMesh.Graphs
using GCMMesh.Partition

FT = Float64
# plotting a rectangular mesh
x1min, x1max = FT(0), FT(1)
x2min, x2max = FT(0), FT(1)
n1, n2 = 3, 4
per = (false, false)

msh1 = equispaced_rectangular_mesh(x1min, x2min, x1max, x2max, n1, n2, per)

scatter(msh1.coordinates[:, 1], msh1.coordinates[:, 2], legend = false)

for i = 1:msh1.nfaces
    plot!(
        msh1.coordinates[msh1.face_verts[i, :], 1],
        msh1.coordinates[msh1.face_verts[i, :], 2],
        legend = false,
    )
end

savefig("rectangular_mesh.png")
# plotting a sphere mesh
ne = 4
radius = FT(10)
msh = sphere_mesh(ne, radius, EquiangularSphereWarp())
#msh = sphere_mesh(ne, FT, EquidistantSphereWarp())
#msh = cube_panel_mesh(ne, FT)


scatter(msh.coordinates[:, 1], msh.coordinates[:, 2], msh.coordinates[:, 3], legend = false)

for i = 1:msh.nfaces
    plot!(
        msh.coordinates[msh.face_verts[i, :], 1],
        msh.coordinates[msh.face_verts[i, :], 2],
        msh.coordinates[msh.face_verts[i, :], 3],
        legend = false,
        aspect_ratio = 1,
    )
end

savefig("sphere_mesh.png")

# plotting partitions for a rectangular mesh
FT = Float64

x1min, x1max = FT(0), FT(1)
x2min, x2max = FT(0), FT(1)
limits = (x1min, x2min, x1max, x2max)
n1, n2 = 16, 16 # # of elements along x1 and x2 axis respectively
#n1, n2 = 4, 4 # # of elements along x1 and x2 axis respectively
#npart = 4 # # of partitions
#npart = 5 # # of partitions
#npart = 2 # # of partitions
npart = 3 # # of partitions

mesh = equispaced_rectangular_mesh(limits..., n1, n2, (false, false))

partition = rsb_partition(mesh, npart)
colors = (:red, :blue, :green, :yellow, :brown)

scatter(mesh.coordinates[:, 1], mesh.coordinates[:, 2], legend = false)

for i = 1:mesh.nelems
    verts = mesh.elem_verts[i, :]
    x = mesh.coordinates[verts, 1]
    y = mesh.coordinates[verts, 2]
    color = colors[partition[i]+1]
    shape = Shape(x, y)
    plot!(shape, c = color)
end

#savefig("rectangle_mesh_4_partitions.png")
#savefig("rectangle_mesh_5_partitions.png")
#savefig("rectangle_mesh_2_partitions.png")
savefig("rectangle_mesh_3_partitions.png")

# plotting a sphere mesh
colors = (:red, :blue, :green, :yellow, :brown)
ne = 4
radius = FT(10)
mesh = cube_panel_mesh(ne, FT)

#npart = 4
npart = 2
#npart = 3
partition = rsb_partition(mesh, npart)

elems_per_panel = div(mesh.nelems, 6)

plot()

for i = 1:mesh.nelems
    verts = mesh.elem_verts[i, :]
    x = mesh.coordinates[verts, 1]
    y = mesh.coordinates[verts, 2]
    z = mesh.coordinates[verts, 3]
    color = colors[partition[i]+1]
    panelno = cld(i, elems_per_panel)
    xx, yy = unfold_cube_panel_to_plane(x, y, z, panelno)
    plot!(Shape(xx, yy), c = color, legend = false)
end

savefig("sphere_mesh_2_partitions.png")
#savefig("sphere_mesh_3_partitions.png")
#savefig("sphere_mesh_4_partitions.png")
