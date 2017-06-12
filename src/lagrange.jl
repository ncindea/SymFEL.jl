"""
    get_lagrange_basis(n = 1, varcoeff = false)

Get Lagrange basis function of order `n`.
"""
function get_lagrange_basis(n = 1, varcoeff = false)
  E = eye(Int64, n+1, n+1)
  symString = ""
  for i = 1:n
    symString *= ("a" * bin(i - 1) * " ")
  end
  symString *= ("a" * bin(n))

  L = symbols(symString)

  basis = []
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
    Leq = []

    for j = 1:(n + 1)
      eq = subs(p, x, points[j])
      push!(Leq, eq)
    end

    Leq = Leq - E[1:(n+1), i]
    LL = solve(convert(Array{SymPy.Sym, 1}, Leq), collect(L[1:(n+1)]))

    q = 0
    for j = 1:(n+1)
      q = q + LL[L[j]] * x^(j - 1)
    end
    push!(basis, q)
  end
  basis
end

"""
    get_lagrange_em(p = 1, m = 0, n = 0)

Get Lagrange finite elements elementary matrices.

# Arguments
  * `p`: degree of polynomials in the basis.
  * `m`: number of derivatives on the first function.
  * `n`: number of derivatives on the second function.
"""
function get_lagange_em(p = 1, m = 0, n = 0)
  M = Array{SymPy.Sym}(p+1, p+1)
  F = get_lagrange_basis(p)
  for i = 1:p+1
    for j = 1:p+1
      M[i, j] = simplify(integrate(diff(F[i], x, m) * diff(F[j], x, n), (x, 0, h)))
    end
  end
  M
end

"""
    get_lagrange_em_varcoeff(p = 1, m = 0, n = 0, f = 1)
    Get Hermite finite elements elementary matrices for variable coefficients.

# Arguments
  * `p`: degree of polynomials in the basis.
  * `m`: number of derivatives on the first function.
  * `n`: number of derivatives on the second function.
  * `f`: the variable coefficient.
"""
function get_lagrange_em_varcoeff(p = 1, m = 0, n = 0, f = 1)
  M = Array{SymPy.Sym}(p+1, p+1)
  F = get_lagrange_basis(p, true)
  for i = 1:p+1
    for j = 1:p+1
      M[i, j] = simplify(integrate(diff(F[i], x, m) * diff(F[j], x, n) * f, (x, xa, xb)))
      M[i, j] = simplify(subs(M[i, j], h, xb - xa))
    end
  end
  M
end
