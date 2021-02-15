# FEMTools.jl Documentation

This package contains several functions usefull for the implementation
of the finite elements method (FEM). More preciselly, the package
`SymPy` is used
in order to compute exactly the elementary finite element functions for both
Lagrange and Hermite finite elements.

```@docs
FEMTools.get_em(deg1=1, deg2=1, der1=0, der2=0; fe1="Lagrange", fe2="Lagrange")
```

## Lagrange finite elements

```@docs
FEMTools.get_lagrange_basis(n = 1, varcoeff = false)
```

```@docs
FEMTools.get_lagrange_em(p = 1, m = 0, n = 0)
```

```@docs
FEMTools.get_lagrange_em_varcoeff(p = 1, m = 0, n = 0, f = 1)
```

```@docs
FEMTools.get_square_lagrange_em((px, py) = (1, 1), (mx, my) = (0, 0), (nx, ny) = (0, 0))
```
    


## Hermite finite elements

```@docs
FEMTools.get_hermite_basis(n = 3, varcoeff = false)
```

```@docs
FEMTools.get_hermite_em(p = 3, m = 0, n = 0)
```

```@docs
FEMTools.get_hermite_em_varcoeff(p = 3, m = 0, n = 0, f = 1)
```

```@docs
FEMTools.interpolate(fd, t, ti)
```

## Assembling functions -- 1d

```@docs
FEMTools.assemble_1d_FE_matrix(elem::Array{Float64, 2}, nbNodes::Int64;
    intNodes = 0, dof1 = 1, dof2 = 1)
```

```@docs
FEMTools.assemble_1d_nu_FE_matrix(elem::Matrix{SymPy.Sym}, nodes::Array{Float64, 1};
    intNodes = 0, dof1 = 1, dof2 = 1)
```

## Assembling functions -- 2d

```@docs
FEMTools.assemble_squaremesh_FE_matrix(el_mat::Array{Float64, 2},
                                       elements::Array{Int64, 2};
                                       order1 = 1,
                                       order2 = 1,
                                       dof1 = 1,
                                       dof2 = 1)
```
