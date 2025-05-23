using SymFEL
using Test
using SymPy
using LinearAlgebra

println("# Testing assemble multiplicative coeff in 1d")
elet = SymFEL.get_etensor(1, 1, 1, 0, 0, 0)

nodes = convert(Array{Float64, 1}, range(0, stop=1, length=101))
dx = nodes[2] - nodes[1]

elet_dx = zeros(size(elet))
for i = 1:size(elet, 3)
    elet_dx[:,:,i] = elet[:,:,i].subs(h, dx)
end


a = cos.(nodes)
M = SymFEL.assemble_1d_FE_matrix_multcoeff(elet_dx,
                                           a,
                                           101)

x = ones(101)
y = sin.(nodes)
@test abs(x' * M * y - (1 - cos(2))/4) <= 1e-4
