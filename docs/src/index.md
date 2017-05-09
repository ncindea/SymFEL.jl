# FEMTools Documentation

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
FEMTools.get_lagange_em(p = 1, m = 0, n = 0)
```

```@docs
FEMTools.get_lagrange_em_varcoeff(p = 1, m = 0, n = 0, f = 1)
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

## Assembling functions

```@docs
FEMTools.assemble_1d_FE_matrix(elem, nbNodes)
```
