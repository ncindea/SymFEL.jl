using FE
using Test
using SymPy
using LinearAlgebra
x = SymPy.symbols("x")
h = SymPy.symbols("h")
xa, xb = SymPy.symbols("xa xb")

# Testing Lagrange elements
println("# Testing Lagrange elements")
p1 = FE.get_lagrange_basis()
@test SymPy.degree(p1[1], gen=x) == 1
@test SymPy.degree(p1[2], gen=x) == 1
@test SymPy.subs(p1[1], x, 0) == 1
@test SymPy.subs(p1[1], x, h) == 0
@test SymPy.subs(p1[2], x, h) == 1
@test SymPy.subs(p1[2], x, 0) == 0

p4 = FE.get_lagrange_basis(4)
E = Matrix(1I, 5, 5)
for i = 1:4
  @test SymPy.degree(p4[i], gen=x) == 4
  for j = 1:5
      @test SymPy.subs(p4[i], x, (j - 1) * h / 4) == E[i, j]
  end
end

M = FE.get_lagrange_em()
@test (M[1, 1] == h / 3) && ((M[1, 2] == h / 6)) &&
      (M[2, 1] == h / 6) && ((M[2, 2] == h / 3))

MVC = FE.get_lagrange_em_varcoeff()
MVC = MVC.subs([(xa, 0), (xb, h)])
@test (MVC[1, 1] == h / 3) && (MVC[1, 2] == h / 6) &&
      (MVC[2, 1] == h / 6) && (MVC[2, 2] == h / 3)

# Testing Hermite elements
println("# Testing Hermite elements")
p3 = FE.get_hermite_basis()
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

MH = FE.get_hermite_em()
MH_hc = [13*h/35 11*h^2/210 9*h/70 -13*h^2/420;
         11*h^2/210 h^3/105 13*h^2/420 -h^3/140;
         9*h/70 13*h^2/420 13*h/35 -11*h^2/210;
         -13*h^2/420 -h^3/140 -11*h^2/210 h^3/105]
@test MH == MH_hc

MHVC = FE.get_hermite_em_varcoeff()
MHVC = MHVC.subs([(xa, 0), (xb, h)])
@test MHVC == MH_hc

try
  FE.get_hermite_basis(2)
catch y
  println(y)
end

##
t = range(0, stop=1,length=11)
t = convert(Array{Float64}, t)
ti = range(0, stop=1, length=101)
x = sin.(pi * t')
y = pi * cos.(pi * t')
fdf = [x; y]
f = FE.interpolate(fdf, t, ti)
@test norm(f[1:10:101]' - x) <= 1e-15

## Testing get_em() function.
println("Testing get_em() function.")
@test FE.get_em(1, 1, 0, 0) == FE.get_lagrange_em(1, 0, 0)
@test FE.get_em(3, 3, 1, 0; fe1 = "Hermite", fe2 = "Hermite") == FE.get_hermite_em(3, 1, 0)
try
  FE.get_em(1, 1, 0, 0; fe1 = "a", fe2 = "Lagrange")
catch y
  println(y)
end

try
  FE.get_em(1, 1, 2, 0)
catch y
  println(y)
end

try
  FE.get_em(1, 1, 1, 2)
catch y
  println(y)
end

try
  FE.get_em(3, 3, 3, 0; fe1 = "Hermite", fe2="Hermite")
catch y
  println(y)
end

try
  FE.get_em(3, 3, 2, -1; fe1 = "Hermite", fe2="Hermite")
catch y
  println(y)
end

##
n = 1001
x = range(0, stop=1, length=2*n-1).^2
dx = 1 / (n - 1)

elem_M = FE.get_lagrange_em(2, 0, 0)
elem_M = elem_M.subs(h, dx)
elem_M = convert(Matrix{Float64}, elem_M)
M = FE.assemble_1d_FE_matrix(elem_M, n, intNodes = 1, dof1 = 1, dof2 = 1)
@test abs(convert(Float64, (x' * M * x)[1]) - 1 / 5) < 1e-15

elem_M_p1 = FE.get_lagrange_em(1, 0, 0)
elem_M_p1 = elem_M_p1.subs(h, dx)
elem_M_p1 = convert(Matrix{Float64}, elem_M_p1)
M_p1 = FE.assemble_1d_FE_matrix(elem_M_p1, n)
x = range(0, stop=1, length=n).^2
@test abs(convert(Float64, (x' * M_p1 * x)[1]) - 1 / 5) < 1e-6

elem_M_h3 = FE.get_hermite_em(3, 0, 0)
elem_M_h3 = elem_M_h3.subs(h, dx)
elem_M_h3 = convert(Matrix{Float64}, elem_M_h3)
M_h3 = FE.assemble_1d_FE_matrix(elem_M_h3, n, intNodes = 0, dof1 = 2, dof2 = 2)
x = zeros(2 * n, 1)
x[1:2:(2*(n-1)+1)] = convert(Array{Float64}, range(0, stop=1, length=n).^3)
x[2:2:(2*(n-1)+2)] = 3*convert(Array{Float64}, range(0, stop=1, length=n).^2)
@test abs(convert(Float64, (x' * M_h3 * x)[1]) - 1 / 7) < 1e-15

nodes = convert(Array{Float64, 1}, range(0, stop=1, length=n))
elem_M_p1 = FE.get_lagrange_em(1, 0, 0)
M_p1_nu = FE.assemble_1d_nu_FE_matrix(elem_M_p1, nodes; intNodes = 0, dof1 = 1, dof2 = 1)
A = M_p1 - M_p1_nu
@test norm(A[:]) < 1e-15
##

nodes = [convert(Array{Float64}, range(0, stop=0.99, length=100));
 convert(Array{Float64}, range(0.99, stop=1, length=101))]
M_p1_nu2 = FE.assemble_1d_nu_FE_matrix(elem_M_p1, nodes; intNodes = 0, dof1 = 1, dof2 = 1)
x = nodes.^2
@test abs((x' * M_p1_nu2 * x)[1] - 1 / 5) < 1e-4
println("All tests passed.")
