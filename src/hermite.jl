"""
    get_hermite_basis(n = 3, varcoeff = false)

Get Hermite basis function of order `n`.
"""
function get_hermite_basis(n = 3, varcoeff = false)
  if (n < 3) || (n % 2 != 1)
    error("The degree of Hermite polynomials should be odd and >= 3")
  end
  E = eye(Int64, n + 1, n + 1)
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

  P = [] # the list containing the derivatives of p
  push!(P, p)
  for i = 1:((n+1)/2 - 1)
    push!(P, diff(p, x, convert(Int64, i)))
  end

  for i = 1:(n + 1)
    Leq = [];
    for j = 1:convert(Int64, (n+1)/2 )
      if (varcoeff)
        push!(Leq, subs(P[j], x, xa))
      else
        push!(Leq, subs(P[j], x, 0))
      end
    end
    for j = 1:convert(Int64, (n+1)/2 )
      if (varcoeff)
        push!(Leq, subs(P[j], x, xb))
      else
        push!(Leq, subs(P[j], x, h))
      end
    end
    Leq = Leq - E[1:(n+1), i]
    LL = solve(convert(Array{SymPy.Sym, 1}, Leq), collect(L[1:(n+1)]))

    p = 0
    for j = 1:(n+1)
      p = p + LL[L[j]] * x^(j - 1)
    end
    push!(basis, p)
  end
  basis
end


"""
    get_hermite_em(p = 3, m = 0, n = 0)

Get Hermite finite elements elementary matrices.

# Arguments
  * `p`: degree of polynomials in the basis.
  * `m`: number of derivatives on the first function.
  * `n`: number of derivatives on the second function.
"""
function get_hermite_em(p = 3, m = 0, n = 0)
  M = Array{SymPy.Sym}(p+1, p+1)
  F = get_hermite_basis(p)
  for i = 1:p+1
    for j = 1:p+1
      M[i, j] = simplify(integrate(diff(F[i], x, m) * diff(F[j], x, n), (x, 0, h)))
    end
  end
  M
end

"""
    get_hermite_em_varcoeff(p = 3, m = 0, n = 0, f = 1)
    Get Hermite finite elements elementary matrices for variable coefficients.

# Arguments
  * `p`: degree of polynomials in the basis.
  * `m`: number of derivatives on the first function.
  * `n`: number of derivatives on the second function.
  * `f`: the variable coefficient.
"""
function get_hermite_em_varcoeff(p = 3, m = 0, n = 0, f = 1)
  M = Array{SymPy.Sym}(p+1, p+1)
  F = get_hermite_basis(p, true)
  for i = 1:p+1
    for j = 1:p+1
      M[i, j] = simplify(integrate(diff(F[i], x, m) * diff(F[j], x, n) * f, (x, xa, xb)))
    end
  end
  M
end

"""
    interpolate(fd::Matrix{Float64}, t::Vector{Float64}, ti::Vector{Float64})

Interpolates `fd` from `t` to `ti`.

# Arguments:
  * `fd`: values and derivatives of function to interpolate.
  * `t`: values in which `fd` is known.
  * `ti`: values in which `fd` is interpolated.
"""
function interpolate(fd, t, ti)
  p = size(fd, 1) # 2*p -1 is the degree of polynomials
  n = size(fd, 2) # number of samples
  if (n != length(t))
    error("The number of columns in the Matrix fd should be egal to the number of elements in the Vector t.")
  end
  if (norm(diff(diff(t))) > 1e-12 || norm(diff(diff(ti))) > 1e-12 ||
        abs(t[1]-ti[1])>1e-12 || abs(t[end] - ti[end]) > 1e-12)
    error("Vectors t and ti should be uniform partitions of [a,b] interval.")
  end
  F = get_hermite_basis(2*p - 1)
  DT = t[2] - t[1];
  FDT = []
  for i = 1:2*p
    push!(FDT, SymPy.subs(F[i], h, DT))
  end
  ni = length(ti)
  fi = collect(zeros(ni , 1))
  for i = 1:ni-1
    j = floor(ti[i] / DT) + 1
    j = convert(Int64, j)
    xx = ti[i] - t[j]
    S = 0;
    for l = 1:p
      S = S + fd[l, j] * SymPy.subs(FDT[l], x, xx)
      S = S + fd[l, j + 1] * SymPy.subs(FDT[l+p], x, xx)
    end
    fi[i] = S
    if (isnan(S))
      fi[i] = 0
    end

  end
  fi[ni] = fd[1, n]
  fi
end
