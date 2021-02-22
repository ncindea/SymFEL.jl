# Copyright 2017-2021 Nicolae Cindea
# See accompanying license file.

module SymFEL
__precompile__(false)

using SymPy

# polynomial variable

include("lagrange.jl")
include("hermite.jl")
include("assemble1d.jl")
include("assemble2d.jl")


export x
export y
export z
export h
export xa, xb

"""
    x = SymPy.symbols("x")

Symbolic variable x.

This variable is exported.
"""
x = SymPy.symbols("x")

"""
    y = SymPy.symbols("y")

Symbolic variable y.

This variable is exported.
"""
y = SymPy.symbols("y")

"""
    z = SymPy.symbols("z")

Symbolic variable z.

This variable is exported.
"""
z = SymPy.symbols("z")

"""
    h = SymPy.symbols("h")

Symbolic variable h.

This variable is exported.
"""
h = SymPy.symbols("h")

"""
    xa, xb = SymPy.symbols("xa xb")

Symbolic vaiables xa and xb.

This variables are used for variable coefficients problems. This variables are exported.
"""
xa, xb = SymPy.symbols("xa xb")

"""
    get_em(deg1=1, deg2=1, der1=0, der2=0; fe1="Lagrange", fe2="Lagrange")

Gives 1d elementary matrices for Lagrange and Hermite elements.

The elementary matrices are computed for elements of length `h`.
"""
function get_em(deg1=1, deg2=1, der1=0, der2=0; fe1="Lagrange", fe2="Lagrange")
    global x
    global h
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
    M = Array{SymPy.Sym}(undef, l1, l2)
    for i = 1:l1
        for j = 1:l2
            M[i, j] = simplify(integrate(diff(p1[i], x, der1) * diff(p2[j], x, der2), (x, 0, h)))
        end
    end
    M
end


end # module
