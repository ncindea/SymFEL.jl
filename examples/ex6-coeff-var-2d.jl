using Revise
using SymPy
using SymFEL
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
# the file is square.jld2 is preparedseparetely
k = @load "square-simple.jld2"

x, h = symbols("x h")

elem_M = SymFEL.get_square_lagrange_em((1, 1), (0, 0), (0, 0))

elet_M = SymFEL.get_etensor(1, 1, 1, 0, 0, 0)
elet_K = SymFEL.get_etensor(1, 1, 1, 1, 1, 0)

nc = ([1, 2, 2, 1],
      [1, 1, 2, 2])
el_ten_K = (SymFEL.get_square_etensor(elet_K,
                                      elet_M,
                                      nc,
                                      nc,
                                      nc)
            + SymFEL.get_square_etensor(elet_M,
                                      elet_K,
                                      nc,
                                      nc,
                                      nc))

dx = norm(nodes_coordinate[:, elements_bound[1, 1]] -
    nodes_coordinate[:, elements_bound[2, 1]])
elem_M_dx = convert(Matrix{Float64}, elem_M.subs(h, dx))
el_ten_K_dx = zeros(size(el_ten_K))
for k = 1:size(el_ten_K, 3)
    el_ten_K_dx[:, :, k] = el_ten_K[:, :, k].subs(h, dx)
end


M = SymFEL.assemble_squaremesh_FE_matrix(elem_M_dx, elements_int)

xx = nodes_coordinate[1,:]
yy = nodes_coordinate[2,:]
a = 1 .+ xx .* yy

sx = sin.(pi * xx)
sy = sin.(pi * yy)
cx = cos.(pi * xx)
cy = cos.(pi * yy)

F = @. (-pi * sy * (yy * cx - pi * a * sx)
        -pi * sx * (xx * cy - pi * a * sy)
        + sx * sy)

uex = sx .* sy

K = SymFEL.assemble_squaremesh_FE_matrix_coeffmult(el_ten_K_dx, a, elements_int)


A = K + M

# boundary condition
tgv = 1e30
A[nodes_boundary, nodes_boundary] += tgv * sparse(Matrix{Float64}(I, nodes_boundary_N, nodes_boundary_N))

u = A \ (M * F)

eru = u - uex
er = sqrt(eru' * M * eru)
println(er)
