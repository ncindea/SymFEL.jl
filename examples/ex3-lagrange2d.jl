# # An equation of second order in two dimmensions
#
# Let \\(\Omega = (0, 1)^2\\) the unit square and denote \\( \Gamma = \partial \Omega\\) its boundary.
# Consider the following problem. Given \\(f \in C(\Omega)\\), find a function \\(u\\) satisfying
# -\Delta u + u = f in \Omega
# u = 0 on \Gamma

using SymFEL
using SymPy
using LinearAlgebra
using SparseArrays
using WriteVTK
using JLD2

## discretization parameters
# we use the mesh square.msh (in gmsh format)
# obtained from square.geo using gmsh
# this mesh is formed by quad elements
# we use gmsh module for read the mesh
# there are some problems using Threads and gmsh (in linux)
# the file is square.jld2 is prepared separetely
@load "square.jld2"

# elementary matrices - P2 x P2
elem_Mxy = SymFEL.get_square_lagrange_em((2, 2), (0, 0), (0, 0))
elem_Kxy = SymFEL.get_square_lagrange_em((2, 2), (1, 0), (1, 0)) +
    SymFEL.get_square_lagrange_em((2, 2), (0, 1), (0, 1))

dx = norm(nodes_coordinate[:, elements_bound[1,1]] - nodes_coordinate[:, elements_bound[2,1]])

elem_Kxy_dx = convert(Matrix{Float64}, elem_Kxy.subs(h, dx))
elem_Mxy_dx = convert(Matrix{Float64}, elem_Mxy.subs(h, dx));

# global matrices
K = SymFEL.assemble_squaremesh_FE_matrix(elem_Kxy_dx, elements_int, order1=2, order2=2)
M = SymFEL.assemble_squaremesh_FE_matrix(elem_Mxy_dx, elements_int, order1=2, order2=2)

f = (2*pi^2 + 1) * sin.(pi * nodes_coordinate[1,:]) .* sin.(pi * nodes_coordinate[2,:])

F = M * f

A = K + M

# boundary condition
tgv = 1e30
A[nodes_boundary, nodes_boundary] += tgv * sparse(Matrix{Float64}(I, nodes_boundary_N, nodes_boundary_N))

u = A \ F
u_exact = sin.(pi*nodes_coordinate[1,:]) .* sin.(pi*nodes_coordinate[2,:])
err = u - u_exact

println("L2 error : ", sqrt(err' * M * err))
println("H1 error : ", sqrt(err' * K * err))


# export to vtk
cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_QUAD, elements_int[1:8, i]) for i = 1:elements_int_N]


points_x = nodes_coordinate[1, :]
points_y = nodes_coordinate[2, :]
vtkfile = vtk_grid("ex3-output", points_x, points_y, cells)

vtkfile["my_point_data", VTKPointData()] = u
outfiles = vtk_save(vtkfile)
