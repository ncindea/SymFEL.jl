# # An equation of second order in two dimmensions
#
# Let \\(\Omega = (0, 1)^2\\) the unit square and denote \\( \Gamma = \partial \Omega\\) its boundary.
# Consider the following problem. Given \\(f \in C(\Omega)\\), find a function \\(u\\) satisfying
# -\Delta u + u = f in \Omega
# u = 0 on \Gamma

using Revise
using FEMTools
using SymPy
using LinearAlgebra
using SparseArrays
import gmsh
using WriteVTK


#using PyPlot
#close("all")

# symbols 
x = SymPy.symbols("x")
h = SymPy.symbols("h")

## discretization parameters
# we use the mesh square.msh (in gmsh format)
# obtained from square.geo using gmsh
# this mesh is formed by quad elements
# we use gmsh module for read the mesh
gmsh.initialize()
gmsh.open("square.msh")
nodes = gmsh.model.mesh.getNodes()
nodes_label = convert(Array{Int64}, nodes[1])
nodes_N = length(nodes_label)
nodes_coordinate = reshape(nodes[2], (3, nodes_N))
nodes_boundary = convert(Array{Int64}, unique(gmsh.model.mesh.getNodesByElementType(1)[1]))
nodes_boundary_N = length(nodes_boundary)

elements = gmsh.model.mesh.getElements()

elements_bound_label = convert(Array{Int64}, elements[2][1])
elements_bound_N = length(elements_bound_label)

elements_int_label = convert(Array{Int64}, elements[2][2])
elements_int_N = length(elements_int_label)

elements_bound = reshape(convert(Array{Int64}, elements[3][1]), (3, elements_bound_N))
elements_int = reshape(convert(Array{Int64}, elements[3][2]), (9, elements_int_N))
#gmsh.fltk.run()
gmsh.finalize()

# elementary matrices - P2 x P2
elem_Mxy = FEMTools.get_square_lagrange_em((2, 2), (0, 0), (0, 0))
elem_Kxy = FEMTools.get_square_lagrange_em((2, 2), (1, 0), (1, 0)) +
    FEMTools.get_square_lagrange_em((2, 2), (0, 1), (0, 1))

dx = norm(nodes_coordinate[:, elements_bound[1,1]] - nodes_coordinate[:, elements_bound[2,1]])

elem_Kxy_dx = convert(Matrix{Float64}, elem_Kxy.subs(h, dx))
elem_Mxy_dx = convert(Matrix{Float64}, elem_Mxy.subs(h, dx));

# global matrices
K = FEMTools.assemble_squaremesh_FE_matrix(elem_Kxy_dx, elements_int, order1=2, order2=2)
M = FEMTools.assemble_squaremesh_FE_matrix(elem_Mxy_dx, elements_int, order1=2, order2=2)

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
