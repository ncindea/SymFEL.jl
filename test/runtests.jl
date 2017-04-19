using FEMTools
using Base.Test

# Testing Lagrange elements
println("# Testing Lagrange elements")
p1 = FEMTools.get_lagrange_basis()
@test SymPy.degree(p1[1], gen=FEMTools.x) == 1
@test SymPy.degree(p1[2], gen=FEMTools.x) == 1
@test SymPy.subs(p1[1], FEMTools.x, 0) == 1
@test SymPy.subs(p1[1], FEMTools.x, FEMTools.h) == 0
@test SymPy.subs(p1[2], FEMTools.x, FEMTools.h) == 1
@test SymPy.subs(p1[2], FEMTools.x, 0) == 0

p4 = FEMTools.get_lagrange_basis(4)
E = eye(5, 5)
for i = 1:4
  @test SymPy.degree(p4[i], gen=FEMTools.x) == 4
  for j = 1:5
      @test SymPy.subs(p4[i], FEMTools.x, (j - 1) * FEMTools.h / 4) == E[i, j]
  end
end

M = FEMTools.get_lagange_em()
@test (M[1, 1] == FEMTools.h / 3) && ((M[1, 2] == FEMTools.h / 6)) &&
      (M[2, 1] == FEMTools.h / 6) && ((M[2, 2] == FEMTools.h / 3))

MVC = FEMTools.get_lagrange_em_varcoeff()
MVC = SymPy.subs(MVC, FEMTools.xa, 0)
MVC = SymPy.subs(MVC, FEMTools.xb, FEMTools.h)
@test (MVC[1, 1] == FEMTools.h / 3) && (MVC[1, 2] == FEMTools.h / 6) &&
      (MVC[2, 1] == FEMTools.h / 6) && (MVC[2, 2] == FEMTools.h / 3)

# Testing Hermite elements
println("# Testing Hermite elements")
p3 = FEMTools.get_hermite_basis()
@test SymPy.subs(p3[1], FEMTools.x, 0) == 1
@test SymPy.subs(p3[1], FEMTools.x, FEMTools.h) == 0
@test SymPy.subs(SymPy.diff(p3[1], FEMTools.x), FEMTools.x, 0) == 0
@test SymPy.subs(SymPy.diff(p3[1], FEMTools.x), FEMTools.x, FEMTools.h) == 0

@test SymPy.subs(p3[2], FEMTools.x, 0) == 0
@test SymPy.subs(p3[2], FEMTools.x, FEMTools.h) == 0
@test SymPy.subs(SymPy.diff(p3[2], FEMTools.x), FEMTools.x, 0) == 1
@test SymPy.subs(SymPy.diff(p3[2], FEMTools.x), FEMTools.x, FEMTools.h) == 0

@test SymPy.subs(p3[3], FEMTools.x, 0) == 0
@test SymPy.subs(p3[3], FEMTools.x, FEMTools.h) == 1
@test SymPy.subs(SymPy.diff(p3[3], FEMTools.x), FEMTools.x, 0) == 0
@test SymPy.subs(SymPy.diff(p3[3], FEMTools.x), FEMTools.x, FEMTools.h) == 0

@test SymPy.subs(p3[4], FEMTools.x, 0) == 0
@test SymPy.subs(p3[4], FEMTools.x, FEMTools.h) == 0
@test SymPy.subs(SymPy.diff(p3[4], FEMTools.x), FEMTools.x, 0) == 0
@test SymPy.subs(SymPy.diff(p3[4], FEMTools.x), FEMTools.x, FEMTools.h) == 1

MH = FEMTools.get_hermite_em()
h = FEMTools.h
MH_hc = [13*h/35 11*h^2/210 9*h/70 -13*h^2/420;
         11*h^2/210 h^3/105 13*h^2/420 -h^3/140;
         9*h/70 13*h^2/420 13*h/35 -11*h^2/210;
         -13*h^2/420 -h^3/140 -11*h^2/210 h^3/105]
@test MH == MH_hc

MHVC = FEMTools.get_hermite_em_varcoeff()
MHVC = SymPy.subs(MHVC, FEMTools.xa, 0)
MHVC = SymPy.subs(MHVC, FEMTools.xb, FEMTools.h)
@test MHVC == MH_hc

println("All tests passed.")
