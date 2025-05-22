using LinearAlgebra
"""
        get_lagrange_basis(n = 1, varcoeff = false;
                           x = symbols("x"), h = symbols("h"))

    Get Lagrange basis function of order `n`.
    """
function get_lagrange_basis(n = 1, varcoeff = false;
                            x = symbols("x"), h = symbols("h"))
    E = Matrix(1I, n+1, n+1)
    symString = ""
    for i = 1:n
        symString *= ("a" * string(i - 1, base=10) * " ")
    end
    symString *= ("a" * string(n, base=10))

    L = SymPy.symbols(symString)

    basis = Array{SymPy.Sym}(undef, n+1)
    # p is a polynomial of degree n
    p = 0

    for i = 1:n+1
        p = p + L[i]*x^(i-1)
    end
    points = []
    for i = 0:n
        push!(points, i * h / n)
    end

    
    for i = 1:(n + 1)
        Leq = zeros(Sym, n+1)

        for j = 1:(n + 1)
            eq = SymPy.subs(p, x, points[j])
            Leq[j] = eq
        end

        Leq = Leq .- E[1:(n+1), i]
        if n == 0
            L2 = zeros(Sym, 1)
            L2[1] = L
            L = L2
        end
        LL = SymPy.solve(Leq, collect(L))
        LL = convert(Dict{SymPy.Sym, SymPy.Sym}, LL)
        q = 0
        for j = 1:(n+1)
            q = q + LL[L[j]] * x^(j - 1)
        end
        basis[i] = q
    end
    basis
end

"""
    get_lagrange_em(p = 1, m = 0, n = 0;
                    x = symbols("x"), h = symbols("h"))

Get Lagrange finite elements elementary matrices.

# Arguments
  * `p`: degree of polynomials in the basis.
  * `m`: number of derivatives on the first function.
  * `n`: number of derivatives on the second function.
"""
function get_lagrange_em(p = 1, m = 0, n = 0;
                         x = symbols("x"), h = symbols("h"))
    M = Array{SymPy.Sym}(undef, p+1, p+1)
    F = get_lagrange_basis(p)
    for i = 1:p+1
        for j = 1:p+1
            M[i, j] = simplify(integrate(diff(F[i], x, m) * diff(F[j], x, n), (x, 0, h)))
        end
    end
    M
end

"""
    get_lagrange_em_varcoeff(p = 1, m = 0, n = 0, f = 1;
                             x = symbols("x"), h = symbols("h"),
                             xa = symbols("xa"), xb = symbols(xb))

Get Hermite finite elements elementary matrices for variable coefficients.

# Arguments
* `p`: degree of polynomials in the basis.
* `m`: number of derivatives on the first function.
* `n`: number of derivatives on the second function.
* `f`: the variable coefficient.
"""
function get_lagrange_em_varcoeff(p = 1, m = 0, n = 0, f = 1;
                                  x = symbols("x"), h = symbols("h"),
                                  xa = symbols("xa"), xb = symbols("xb"))
    M = Array{SymPy.Sym}(undef, p+1, p+1)
    ff = subs(f, x, xa + (xb - xa) / h * x)
    F = get_lagrange_basis(p, true)
    for i = 1:p+1
        for j = 1:p+1
            M[i, j] = simplify(integrate(diff(F[i], x, m) * diff(F[j], x, n) * ff * (xb - xa) / h, (x, 0, h)))
            M[i, j] = simplify(subs(M[i, j], h, xb - xa))
        end
    end
    M
end

"""
    get_square_lagrange_em((px, py) = (1, 1), (mx, my) = (0, 0), (nx, ny) = (0, 0);
                           x = symbols("x"), h = symbols("h"))

Get Lagrange finite elements elementary matrices for a squared element.

# Arguments
  * `(px, py)` : degree of polynomials in the basis.
  * `(mx, my)` : number of derivatives on the first function wrt x and y
  * `(nx, ny)` : number of derivatives on the second function wrt x and y

# Remarks
  The current version works only for px, py in {1, 2}
"""
function get_square_lagrange_em((px, py) = (1, 1),
                                (mx, my) = (0, 0),
                                (nx, ny) = (0, 0);
                                x = symbols("x"), h = symbols("h"))
    Mx = get_lagrange_em(px, mx, nx)
    My = get_lagrange_em(py, my, ny)

    node_x = [(1, 2, 2, 1),
              (1, 3, 3, 1, 2, 3, 2, 1, 2)]
    node_y = [(1, 1, 2, 2),
              (1, 1, 3, 3, 1, 2, 3, 2, 2)]


    p = (px + 1)*(py + 1)
    M = Array{SymPy.Sym}(undef, p, p)

    for i = 1:p
        for j = 1:p
            M[i,j] = Mx[node_x[px][i], node_x[px][j]] * My[node_y[py][i], node_y[py][j]]
        end
    end
    M
end


