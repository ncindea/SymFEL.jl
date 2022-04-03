using SymFEL
using Test
using SymPy
using LinearAlgebra


# Testing Lagrange elements
println("# Testing Lagrange elements")
p1 = SymFEL.get_lagrange_basis()
@test SymPy.degree(p1[1], gen=x) == 1
@test SymPy.degree(p1[2], gen=x) == 1
@test SymPy.subs(p1[1], x, 0) == 1
@test SymPy.subs(p1[1], x, h) == 0
@test SymPy.subs(p1[2], x, h) == 1
@test SymPy.subs(p1[2], x, 0) == 0

p4 = SymFEL.get_lagrange_basis(4)
E = Matrix(1I, 5, 5)
for i = 1:4
    @test SymPy.degree(p4[i], gen=x) == 4
    for j = 1:5
        @test SymPy.subs(p4[i], x, (j - 1) * h / 4) == E[i, j]
    end
end

M = SymFEL.get_lagrange_em()
@test (M[1, 1] == h / 3) && ((M[1, 2] == h / 6)) &&
    (M[2, 1] == h / 6) && ((M[2, 2] == h / 3))

f = convert(Sym, 1)
MNU = SymFEL.get_lagrange_em_varcoeff(1, 0, 0, f)
MNU = MNU.subs([(xa, 0), (xb, h)])
@test (MNU[1, 1] == h / 3) && ((MNU[1, 2] == h / 6)) &&
    (MNU[2, 1] == h / 6) && ((MNU[2, 2] == h / 3))

MVC = SymFEL.get_lagrange_em_varcoeff()
MVC = MVC.subs([(xa, 0), (xb, h)])
@test (MVC[1, 1] == h / 3) && (MVC[1, 2] == h / 6) &&
    (MVC[2, 1] == h / 6) && (MVC[2, 2] == h / 3)

# Testing Hermite elements
println("# Testing Hermite elements")
p3 = SymFEL.get_hermite_basis()
@test SymPy.subs(p3[1], x, 0) == 1
@test SymPy.subs(p3[1], x, h) == 0
@test SymPy.subs(SymPy.diff(p3[1], x), x, 0) == 0
@test SymPy.subs(SymPy.diff(p3[1], x), x, h) == 0

@test SymPy.subs(p3[2], x, 0) == 0
@test SymPy.subs(p3[2], x, h) == 0
@test SymPy.subs(SymPy.diff(p3[2], x), x, 0) == 1
@test SymPy.subs(SymPy.diff(p3[2], x), x, h) == 0

@test SymPy.subs(p3[3], x, 0) == 0
@test SymPy.subs(p3[3], x, h) == 1
@test SymPy.subs(SymPy.diff(p3[3], x), x, 0) == 0
@test SymPy.subs(SymPy.diff(p3[3], x), x, h) == 0

@test SymPy.subs(p3[4], x, 0) == 0
@test SymPy.subs(p3[4], x, h) == 0
@test SymPy.subs(SymPy.diff(p3[4], x), x, 0) == 0
@test SymPy.subs(SymPy.diff(p3[4], x), x, h) == 1

MH = SymFEL.get_hermite_em()
MH_hc = [13*h/35 11*h^2/210 9*h/70 -13*h^2/420;
         11*h^2/210 h^3/105 13*h^2/420 -h^3/140;
         9*h/70 13*h^2/420 13*h/35 -11*h^2/210;
         -13*h^2/420 -h^3/140 -11*h^2/210 h^3/105]
@test MH == MH_hc

MHVC = SymFEL.get_hermite_em_varcoeff()
MHVC = MHVC.subs([(xa, 0), (xb, h)])
@test MHVC == MH_hc

try
    SymFEL.get_hermite_basis(2)
catch y
    println(y)
end

##
t = range(0, stop=1,length=11)
t = convert(Array{Float64}, t)
ti = range(0, stop=1, length=101)
xv = sin.(pi * t')
yv = pi * cos.(pi * t')
fdf = [xv; yv]
f = SymFEL.interpolate(fdf, t, ti)
@test norm(f[1:10:101]' - xv) <= 1e-15

## Testing get_em() function.
println("Testing get_em() function.")
@test SymFEL.get_em(1, 1, 0, 0) == SymFEL.get_lagrange_em(1, 0, 0)
@test SymFEL.get_em(3, 3, 1, 0; fe1 = "Hermite", fe2 = "Hermite") == SymFEL.get_hermite_em(3, 1, 0)
try
    SymFEL.get_em(1, 1, 0, 0; fe1 = "a", fe2 = "Lagrange")
catch er
    println(er)
end

try
    SymFEL.get_em(1, 1, 2, 0)
catch er
    println(er)
end

try
    SymFEL.get_em(1, 1, 1, 2)
catch er
    println(er)
end

try
    SymFEL.get_em(3, 3, 3, 0; fe1 = "Hermite", fe2="Hermite")
catch er
    println(er)
end

try
    SymFEL.get_em(3, 3, 2, -1; fe1 = "Hermite", fe2="Hermite")
catch er
    println(er)
end

##
n = 1001
xv = range(0, stop=1, length=2*n-1).^2
dx = 1 / (n - 1)

elem_M = SymFEL.get_lagrange_em(2, 0, 0)
elem_M = elem_M.subs(h, dx)
elem_M = convert(Matrix{Float64}, elem_M)
M = SymFEL.assemble_1d_FE_matrix(elem_M, n, intNodes1 = 1, intNodes2 = 1, dof1 = 1, dof2 = 1)
@test abs(convert(Float64, (xv' * M * xv)[1]) - 1 / 5) < 1e-15

elem_M_p1 = SymFEL.get_lagrange_em(1, 0, 0)
elem_M_p1 = elem_M_p1.subs(h, dx)
elem_M_p1 = convert(Matrix{Float64}, elem_M_p1)
M_p1 = SymFEL.assemble_1d_FE_matrix(elem_M_p1, n)
xv = range(0, stop=1, length=n).^2
@test abs(convert(Float64, (xv' * M_p1 * xv)[1]) - 1 / 5) < 1e-6

elem_M_h3 = SymFEL.get_hermite_em(3, 0, 0)
elem_M_h3 = elem_M_h3.subs(h, dx)
elem_M_h3 = convert(Matrix{Float64}, elem_M_h3)
M_h3 = SymFEL.assemble_1d_FE_matrix(elem_M_h3, n, intNodes1 = 0, intNodes2 = 0, dof1 = 2, dof2 = 2)
xv = zeros(2 * n, 1)
xv[1:2:(2*(n-1)+1)] = convert(Array{Float64}, range(0, stop=1, length=n).^3)
xv[2:2:(2*(n-1)+2)] = 3*convert(Array{Float64}, range(0, stop=1, length=n).^2)
@test abs(convert(Float64, (xv' * M_h3 * xv)[1]) - 1 / 7) < 1e-15

nodes = convert(Array{Float64, 1}, range(0, stop=1, length=n))
elem_M_p1 = SymFEL.get_lagrange_em(1, 0, 0)
M_p1_nu = SymFEL.assemble_1d_nu_FE_matrix(elem_M_p1, nodes; intNodes1 = 0, intNodes2 = 0, dof1 = 1, dof2 = 1)
A = M_p1 - M_p1_nu
@test norm(A[:]) < 1e-15
##

nodes = [convert(Array{Float64}, range(0, stop=0.99, length=100));
         convert(Array{Float64}, range(0.99, stop=1, length=101))]
M_p1_nu2 = SymFEL.assemble_1d_nu_FE_matrix(elem_M_p1, nodes; intNodes1 = 0, intNodes2=0, dof1 = 1, dof2 = 1)
xv = nodes.^2
@test abs((xv' * M_p1_nu2 * xv)[1] - 1 / 5) < 1e-4


# lagrange - hermite
n = 1001
xv = range(0, stop=1, length=2*n-1).^4
dx = 1 / (n - 1)

elem_M = SymFEL.get_em(2, 3, 1, 0, fe1="Lagrange", fe2="Hermite")
elem_M = elem_M.subs(h, dx)
elem_M = convert(Matrix{Float64}, elem_M)
M = SymFEL.assemble_1d_FE_matrix(elem_M, n, intNodes1 = 1, intNodes2 = 0, dof1 = 1, dof2 = 2)

xr = zeros(2*n)
xvr2 = range(0, stop=1, length=n)
i = 1:n
xr[2*(i .- 1) .+ 1] = xvr2
xr[2*(i .- 1) .+ 2] = ones(n)

@test abs((xv' * M * xr)[1] - 4 / 5) < 1e-14

n = 1001
f = SymFEL.x^2
elem_M_x = SymFEL.get_lagrange_em_varcoeff(1, 0, 0, f)
nodes = range(0, stop=1, length=n)
nodes = convert(Vector{Float64}, nodes)
M = SymFEL.assemble_1d_nu_FE_matrix_varcoeff(elem_M_x, nodes)
vec = ones(n)
III = vec' * M * vec

@test abs(III - 1/3) < 1e-2

println("All tests passed.")


