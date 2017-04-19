# Copyright 2017 Nicolae Cindea
# See accompanying license file.
module FEMTools

using SymPy

# polynomial variable
x = symbols("x")
xa = symbols("xa")
xb = symbols("xb")
h = symbols("h")

include("lagrange.jl")
include("hermite.jl")

"""
  get_em(deg1=1, deg2=1, der1=0, der2=0; fe1="Lagrange", fe2="Lagrange")
  Gives an elementary matrices for different elements in each dimension.
"""
function get_em(deg1=1, deg2=1, der1=0, der2=0; fe1="Lagrange", fe2="Lagrange")
  if fe1 == "Lagrange"
    if der1 ∈ collect(0:deg1)
      p1 = get_lagrange_basis(deg1)
    else
      error("The third argument should verify der1 ∈ [0, …, deg1]...")
    end
  end
  if fe2 == "Lagrange"
    if der2 ∈ collect(0:deg2)
      p2 = get_lagrange_basis(deg2)
    else
      error("The fourth parameter should verify der2 ∈ [0, …, deg2]...")
    end
  end
  if fe1 == "Hermite"
    if der1 ∈ collect(0:floor(Int, deg1 / 2 + 1))
      p1 = get_hermite_basis(deg1)
    else
      error("The third parameter should verify der1 ∈ [0, …,int(deg1 / 2 + 1)]...")
    end
  end
  if fe2 == "Hermite"
    if der2 ∈ collect(0:floor(Int, deg2 / 2 + 1))
      p2 = get_hermite_basis(deg2)
    else
      error("The fourth parameter should verify der2 ∈ [0, …,int(deg2 / 2 + 1)]...")
    end
  end
  if (fe1 != "Lagrange" && fe1 != "Hermite") || (fe2 != "Lagrange" && fe2 != "Hermite")
    error("Only Lagrange and Hermite FEM are implemented.")
  end

  l1 = length(p1)
  l2 = length(p2)
  M = Array{SymPy.Sym}(l1, l2)
  for i = 1:l1
    for j = 1:l2
      M[i, j] = simplify(integrate(diff(p1[i], x, der1) * diff(p2[j], x, der2), (x, 0, h)))
    end
  end
  M
end


end # module
