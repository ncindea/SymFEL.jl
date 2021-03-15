# Copyright 2017-2021 Nicolae Cindea
# See accompanying license file.

module SymFEL
#__precompile__(false)

using SymPy

# polynomial variable

include("lagrange.jl")
include("hermite.jl")
include("assemble1d.jl")
include("assemble2d.jl")
include("assemble3d.jl")


# export x
# export y
# export z
# export h
# export hx
# export hy
# export hz
# export xa, xb

# """
#     x = SymPy.symbols("x")

# Symbolic variable x.

# This variable is exported.
# """
# x = SymPy.symbols("x")

# """
#     y = SymPy.symbols("y")

# Symbolic variable y.

# This variable is exported.
# """
# y = SymPy.symbols("y")

# """
#     z = SymPy.symbols("z")

# Symbolic variable z.

# This variable is exported.
# """
# z = SymPy.symbols("z")

# """
#     h = SymPy.symbols("h")

# Symbolic variable h.

# This variable is exported.
# """
# h = SymPy.symbols("h")

# """
#     hx = SymPy.symbols("hx")

# Symbolic variable h.

# This variable is exported.
# """
# hx = SymPy.symbols("hx")

# """
#     hy = SymPy.symbols("hy")

# Symbolic variable hy.

# This variable is exported.
# """
# hy = SymPy.symbols("hy")


# """
#     hz = SymPy.symbols("hz")

# Symbolic variable hz.

# This variable is exported.
# """
# hz = SymPy.symbols("hz")



# """
#     xa, xb = SymPy.symbols("xa xb")

# Symbolic vaiables xa and xb.

# This variables are used for variable coefficients problems. This variables are exported.
# """
# xa, xb = SymPy.symbols("xa xb")

export define_symbols
const define_symbols = :(SymPy.@vars x y z h hx hy hz xa xb)

"""
    get_em(deg1=1, deg2=1, der1=0, der2=0; fe1="Lagrange", fe2="Lagrange")

Gives 1d elementary matrices for Lagrange and Hermite elements.

The elementary matrices are computed for elements of length `h`.
"""
function get_em(deg1=1, deg2=1, der1=0, der2=0; fe1="Lagrange", fe2="Lagrange")
    eval(define_symbols)
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


"""
    get_square_em(Mx::Array{SymPy.Sym, 2},
                  My::Array{SymPy.Sym, 2},
                  nc::Tuple{Array{Int64,1},Array{Int64,1}},
                  nr::Tuple{Array{Int64,1},Array{Int64,1}})

Get elementary matrices for a squared or a rectangular element.

# Arguments
  * `Mx` : elementary matrix for the x variable.
  * `My` : elementary matrix for the y variable
  * `nc` : order of nodes for column
  * `nr` : order of nodes for row

"""
function get_square_em(Mx::Array{SymPy.Sym, 2},
                       My::Array{SymPy.Sym, 2},
                       nc::Tuple{Array{Int64,1},Array{Int64,1}},
                       nr::Tuple{Array{Int64,1},Array{Int64,1}})
    

    pxc = size(Mx, 2)
    pxr = size(Mx, 1)

    pyc = size(My, 2)
    pyr = size(My, 1)

    pc = pxc * pyc
    pr = pxr * pyr
    
    M = Array{SymPy.Sym}(undef, pr, pc)

    for i = 1:pr
        for j = 1:pc
            M[i,j] = Mx[nr[1][i], nc[1][j]] * My[nr[2][i], nc[2][j]]
        end
    end
    M
end

"""
    get_cube_em(Mx::Array{SymPy.Sym, 2},
                My::Array{SymPy.Sym, 2},
                Mz::Array{SymPy.Sym, 2},
                nc::Tuple{Array{Int64,1},Array{Int64,1},Array{Int64,1}},
                nr::Tuple{Array{Int64,1},Array{Int64,1}},Array{Int64,1})

Get elementary matrices for a cube or a rectangular hexahedron element.

# Arguments
  * `Mx` : elementary matrix for the x variable.
  * `My` : elementary matrix for the y variable
  * `Mz` : elementary matrix for the z variable
  * `nc` : order of nodes for column
  * `nr` : order of nodes for row

"""
function get_cube_em(Mx::Array{SymPy.Sym, 2},
                     My::Array{SymPy.Sym, 2},
                     Mz::Array{SymPy.Sym, 2},
                     nc::Tuple{Array{Int64,1},Array{Int64,1},Array{Int64,1}},
                     nr::Tuple{Array{Int64,1},Array{Int64,1},Array{Int64,1}})
    

    pxc = size(Mx, 2)
    pxr = size(Mx, 1)

    pyc = size(My, 2)
    pyr = size(My, 1)

    pzc = size(Mz, 2)
    pzr = size(Mz, 1)

    pc = pxc * pyc * pzc
    pr = pxr * pyr * pzr
    
    M = Array{SymPy.Sym}(undef, pr, pc)

    for i = 1:pr
        for j = 1:pc
            M[i,j] = Mx[nr[1][i], nc[1][j]] * My[nr[2][i], nc[2][j]] * Mz[nr[3][i], nc[3][j]]
        end
    end
    M
end


end # module

