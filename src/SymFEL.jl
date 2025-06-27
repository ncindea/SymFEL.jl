# Copyright 2017-2025 Nicolae Cindea
# See accompanying license file.

module SymFEL


using SymPy

include("lagrange.jl")
include("hermite.jl")
include("assemble1d.jl")
include("assemble2d.jl")
include("assemble3d.jl")


"""
    get_em(deg1=1, deg2=1, der1::Function, der2::Function;
           fe1="Lagrange", fe2="Lagrange",
           x=symbols("x"), h=symbols("h"))

Gives 1d elementary matrices for Lagrange and Hermite elements.
* der1 and der2 are functions having two arguments:
  - a symbolic expression
  - a symbolic variable

The elementary matrices are computed for elements of length `h`.
"""
function get_em(deg1, deg2, der1::Function, der2::Function;
                fe1="Lagrange", fe2="Lagrange",
                x=symbols("x"), h=symbols("h"))
    if fe1 == "Lagrange"
        p1 = get_lagrange_basis(deg1)
    end
    if fe2 == "Lagrange"
        p2 = get_lagrange_basis(deg2)
    end
    if fe1 == "Hermite"
        p1 = get_hermite_basis(deg1)
    end
    if fe2 == "Hermite"
        p2 = get_hermite_basis(deg2)
    end
    if (fe1 != "Lagrange" && fe1 != "Hermite") || (fe2 != "Lagrange" && fe2 != "Hermite")
        error("Only Lagrange and Hermite FEM are implemented.")
    end

    l1 = length(p1)
    l2 = length(p2)
    M = Array{SymPy.Sym}(undef, l1, l2)
    for i = 1:l1
        for j = 1:l2
            M[i, j] = simplify(integrate(der1(p1[i], x) * der2(p2[j], x), (x, 0, h)))
        end
    end
    M
end

"""
    get_em(deg1=1, deg2=1, der1=0, der2=0;
           fe1="Lagrange", fe2="Lagrange",
           x=symbols("x"), h=symbols("h"))

Gives 1d elementary matrices for Lagrange and Hermite elements.

The elementary matrices are computed for elements of length `h`.
"""
function get_em(deg1=1, deg2=1, der1=0, der2=0;
                fe1="Lagrange", fe2="Lagrange",
                x=symbols("x"), h=symbols("h"))
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
    get_etensor(deg1=1, deg2=1, deg=1, der1=0, der2=0, der3=0;
                fe1="Lagrange", fe2="Lagrange", fe3="Lagrange",
                x=symbols("x"), h=symbols("h"))

Gives 1d elementary tensors for Lagrange and Hermite elements.

This is used, in particular, for ∫ a(x) u(x) ϕ(x) dx.

The elementary matrices are computed for elements of length `h`.
"""
function get_etensor(deg1=1, deg2=1, deg3=1, der1=0, der2=0, der3=0;
                     fe1="Lagrange", fe2="Lagrange", fe3="Lagrange",
                     x=symbols("x"), h=symbols("h"))
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
    if fe3 == "Lagrange"
        if der3 ∈ collect(0:deg3)
            p3 = get_lagrange_basis(deg3)
        else
            error("The fourth parameter should verify der3 ∈ [0, …, deg3]...")
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
    if fe3 == "Hermite"
        if der3 ∈ collect(0:floor(Int, deg3 / 2 + 1))
            p3 = get_hermite_basis(deg3)
        else
            error("The fourth parameter should verify der3 ∈ [0, …,int(deg3 / 2 + 1)]...")
        end
    end
    if (fe1 != "Lagrange" && fe1 != "Hermite") || (fe2 != "Lagrange" && fe2 != "Hermite") || (fe3 != "Lagrange" && fe3 != "Hermite")
        error("Only Lagrange and Hermite FEM are implemented.")
    end

    l1 = length(p1)
    l2 = length(p2)
    l3 = length(p3)
    M = Array{SymPy.Sym}(undef, l1, l2, l3)
    for i = 1:l1
        for j = 1:l2
            for k = 1:l3
                M[i, j, k] = simplify(integrate(diff(p1[i], x, der1) * diff(p2[j], x, der2) * diff(p3[k], x, der3), (x, 0, h)))
            end
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

"""
    get_square_etensor(Tx::Array{SymPy.Sym, 3},
                       Ty::Array{SymPy.Sym, 3},
                       nc::Tuple{Array{Int64,1},Array{Int64,1}},
                       nr::Tuple{Array{Int64,1},Array{Int64,1}},
                       np::Tuple{Array{Int64,1},Array{Int64,1}})

Get elementary matrices for a squared or a rectangular element.

# Arguments
  * `Tx` : elementary tensor for the x variable.
  * `Ty` : elementary tensor for the y variable
  * `nc` : order of nodes for column
  * `nr` : order of nodes for row
  * `np` : order of nodes for depth

"""
function get_square_etensor(Tx::Array{SymPy.Sym, 3},
                            Ty::Array{SymPy.Sym, 3},
                            nc::Tuple{Array{Int64,1},Array{Int64,1}},
                            nr::Tuple{Array{Int64,1},Array{Int64,1}},
                            np::Tuple{Array{Int64,1},Array{Int64,1}}
                            )
    
    pxp = size(Tx, 3)
    pxc = size(Tx, 2)
    pxr = size(Tx, 1)

    pyp = size(Ty, 3)
    pyc = size(Ty, 2)
    pyr = size(Ty, 1)

    pp = pxp * pyp
    pc = pxc * pyc
    pr = pxr * pyr
    
    T = Array{SymPy.Sym}(undef, pr, pc, pp)

    for i = 1:pr
        for j = 1:pc
            for k = 1:pp
                T[i, j, k] = Tx[nr[1][i], nc[1][j], np[1][k]] * Ty[nr[2][i], nc[2][j], np[2][k]]
            end
        end
    end
    T
end



end # module

