# FEMTools

[![Build Status](https://travis-ci.org/ncindea/FEMTools.jl.svg?branch=master)](https://travis-ci.org/ncindea/FEMTools.jl)
[![Coverage Status](https://coveralls.io/repos/ncindea/FEMTools.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/ncindea/FEMTools.jl?branch=master)
[![codecov.io](http://codecov.io/github/ncindea/FEMTools.jl/coverage.svg?branch=master)](http://codecov.io/github/ncindea/FEMTools.jl?branch=master)

This package contains several functions usefull for the implementation of the finite elements method (FEM). More preciselly, the package `SymPy` is used in order to compute exactly the elementary finite element functions for both Lagrange and Hermite finite elements.

<a id='FEMTools.get_em' href='#FEMTools.get_em'>#</a>
**`FEMTools.get_em`** &mdash; *Function*.



get_em(deg1=1, deg2=1, der1=0, der2=0; fe1="Lagrange", fe2="Lagrange")

Gives an elementary matrices for different elements in each dimension.


<a target='_blank' href='https://github.com/ncindea/FEMTools.jl/tree/f4d34ce1f90860bfd59a5343d361aa79ef28b867/src/FEMTools.jl#L17-L21' class='documenter-source'>source</a><br>


<a id='Lagrange-finite-elements-1'></a>

## Lagrange finite elements

<a id='FEMTools.get_lagrange_basis' href='#FEMTools.get_lagrange_basis'>#</a>
**`FEMTools.get_lagrange_basis`** &mdash; *Function*.



```
get_lagrange_basis(n = 1, varcoeff = false)
```

Get Lagrange basis function of order `n`.


<a target='_blank' href='https://github.com/ncindea/FEMTools.jl/tree/f4d34ce1f90860bfd59a5343d361aa79ef28b867/src/lagrange.jl#L1-L5' class='documenter-source'>source</a><br>

<a id='FEMTools.get_lagange_em' href='#FEMTools.get_lagange_em'>#</a>
**`FEMTools.get_lagange_em`** &mdash; *Function*.



```
get_lagrange_em(p = 1, m = 0, n = 0)
```

Get Lagrange finite elements elementary matrices.

**Arguments**

  * `p`: degree of polynomials in the basis.
  * `m`: number of derivatives on the first function.
  * `n`: number of derivatives on the second function.


<a target='_blank' href='https://github.com/ncindea/FEMTools.jl/tree/f4d34ce1f90860bfd59a5343d361aa79ef28b867/src/lagrange.jl#L48-L57' class='documenter-source'>source</a><br>

<a id='FEMTools.get_lagrange_em_varcoeff' href='#FEMTools.get_lagrange_em_varcoeff'>#</a>
**`FEMTools.get_lagrange_em_varcoeff`** &mdash; *Function*.



```
get_lagrange_em_varcoeff(p = 1, m = 0, n = 0, f = 1)
Get Hermite finite elements elementary matrices for variable coefficients.
```

**Arguments**

  * `p`: degree of polynomials in the basis.
  * `m`: number of derivatives on the first function.
  * `n`: number of derivatives on the second function.
  * `f`: the variable coefficient.


<a target='_blank' href='https://github.com/ncindea/FEMTools.jl/tree/f4d34ce1f90860bfd59a5343d361aa79ef28b867/src/lagrange.jl#L69-L78' class='documenter-source'>source</a><br>


<a id='Hermite-finite-elements-1'></a>

## Hermite finite elements

<a id='FEMTools.get_hermite_basis' href='#FEMTools.get_hermite_basis'>#</a>
**`FEMTools.get_hermite_basis`** &mdash; *Function*.



```
get_hermite_basis(n = 3, varcoeff = false)
```

Get Hermite basis function of order `n`.


<a target='_blank' href='https://github.com/ncindea/FEMTools.jl/tree/f4d34ce1f90860bfd59a5343d361aa79ef28b867/src/hermite.jl#L1-L5' class='documenter-source'>source</a><br>

<a id='FEMTools.get_hermite_em' href='#FEMTools.get_hermite_em'>#</a>
**`FEMTools.get_hermite_em`** &mdash; *Function*.



```
get_hermite_em(p = 3, m = 0, n = 0)
```

Get Hermite finite elements elementary matrices.

**Arguments**

  * `p`: degree of polynomials in the basis.
  * `m`: number of derivatives on the first function.
  * `n`: number of derivatives on the second function.


<a target='_blank' href='https://github.com/ncindea/FEMTools.jl/tree/f4d34ce1f90860bfd59a5343d361aa79ef28b867/src/hermite.jl#L60-L69' class='documenter-source'>source</a><br>

<a id='FEMTools.get_hermite_em_varcoeff' href='#FEMTools.get_hermite_em_varcoeff'>#</a>
**`FEMTools.get_hermite_em_varcoeff`** &mdash; *Function*.



```
get_hermite_em_varcoeff(p = 3, m = 0, n = 0, f = 1)
Get Hermite finite elements elementary matrices for variable coefficients.
```

**Arguments**

  * `p`: degree of polynomials in the basis.
  * `m`: number of derivatives on the first function.
  * `n`: number of derivatives on the second function.
  * `f`: the variable coefficient.


<a target='_blank' href='https://github.com/ncindea/FEMTools.jl/tree/f4d34ce1f90860bfd59a5343d361aa79ef28b867/src/hermite.jl#L81-L90' class='documenter-source'>source</a><br>

<a id='FEMTools.interpolate-Tuple{Any,Any,Any}' href='#FEMTools.interpolate-Tuple{Any,Any,Any}'>#</a>
**`FEMTools.interpolate`** &mdash; *Method*.



```
interpolate(fd::Matrix{Float64}, t::Vector{Float64}, ti::Vector{Float64})
```

Interpolates `fd` from `t` to `ti`.

**Arguments:**

  * `fd`: values and derivatives of function to interpolate.
  * `t`: values in which `fd` is known.
  * `ti`: values in which `fd` is interpolated.


<a target='_blank' href='https://github.com/ncindea/FEMTools.jl/tree/f4d34ce1f90860bfd59a5343d361aa79ef28b867/src/hermite.jl#L102-L111' class='documenter-source'>source</a><br>


<a id='Assembling-functions-1'></a>

## Assembling functions

<a id='FEMTools.assemble_1d_FE_matrix-Tuple{Any,Any}' href='#FEMTools.assemble_1d_FE_matrix-Tuple{Any,Any}'>#</a>
**`FEMTools.assemble_1d_FE_matrix`** &mdash; *Method*.



```
assemble_1d_FE_matrix(elem, nbNodes)
```


<a target='_blank' href='https://github.com/ncindea/FEMTools.jl/tree/f4d34ce1f90860bfd59a5343d361aa79ef28b867/src/assemble1d.jl#L2-L4' class='documenter-source'>source</a><br>
