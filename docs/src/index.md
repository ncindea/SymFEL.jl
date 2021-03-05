# SymFEL.jl Documentation



| Testing                                                                                                       |
|---------------------------------------------------------------------------------------------------------------|
| [![Build Status](https://travis-ci.org/ncindea/SymFEL.jl.svg?branch=master)](https://travis-ci.org/ncindea/SymFEL.jl) |


This package contains several functions usefull for the implementation of the finite elements method (FEM). We use the package `SymPy` for the computation of the basis finite element functions for both Lagrange and Hermite finite elements. We also compute the elementary matrices in 1d and for quadragular elements in 2d. 
For the two dimensional meshes we employ [Gmsh](https://gmsh.info/). The Gmsh SDK should be installed in your path following the instruction on the Gmsh site.

The utilisation of the package is illustrated by several examples in 1d and 2d.

A documentation website is available at [http://lmbp.uca.fr/~cindea/software/SymFEL.jl/](http://lmbp.uca.fr/~cindea/software/SymFEL.jl/).

Here are some features to develop in the next versions:
- 3d finite elements on quadragular meshes
- 2d triangular finite elements
- 3d tetrahedral finite elements
- combine the symbolic integration to numerical integration to speed-up the execution.

## Contents
```@contents
Pages = ["index.md", "examples1.md", "examples2.md", "examples3.md"]
```

## SymFEl.jl

```@docs
x
y
z
h
hx
hy
hz
```

```@docs
SymFEL.get_em(deg1=1, deg2=1, der1=0, der2=0; fe1="Lagrange", fe2="Lagrange")
```

```@docs
SymFEL.get_square_em(Mx::Array{SymPy.Sym, 2},
                     My::Array{SymPy.Sym, 2},
                     nc::Tuple{Array{Int64,1},Array{Int64,1}},
                     nr::Tuple{Array{Int64,1},Array{Int64,1}})
```

```@docs
SymFEL.get_cube_em(Mx::Array{SymPy.Sym, 2},
                   My::Array{SymPy.Sym, 2},
                   Mz::Array{SymPy.Sym, 2},
                   nc::Tuple{Array{Int64,1},Array{Int64,1},Array{Int64,1}},
                   nr::Tuple{Array{Int64,1},Array{Int64,1},Array{Int64,1}})
```

## Lagrange finite elements

```@docs
SymFEL.get_lagrange_basis(n = 1, varcoeff = false)
```

```@docs
SymFEL.get_lagrange_em(p = 1, m = 0, n = 0)
```

```@docs
SymFEL.get_lagrange_em_varcoeff(p = 1, m = 0, n = 0, f = 1)
```

```@docs
SymFEL.get_square_lagrange_em((px, py) = (1, 1), (mx, my) = (0, 0), (nx, ny) = (0, 0))
```
    


## Hermite finite elements

```@docs
SymFEL.get_hermite_basis(n = 3, varcoeff = false)
```

```@docs
SymFEL.get_hermite_em(p = 3, m = 0, n = 0)
```

```@docs
SymFEL.get_hermite_em_varcoeff(p = 3, m = 0, n = 0, f = 1)
```

```@doc
SymFEL.get_square_hermite_em((px, py) = (3, 3), (mx, my) = (0, 0), (nx, ny) = (0, 0))
```

```@docs
SymFEL.interpolate(fd, t, ti)
```

## Assembling functions -- 1d

```@docs
SymFEL.Mesh1d
```

```@docs
SymFEL.assemble_1d_FE_matrix(elem::Array{Float64, 2}, nbNodes::Int64;
    intNodes1 = 0, intNodes2 = 0, dof1 = 1, dof2 = 1)
```

```@docs
SymFEL.assemble_1d_nu_FE_matrix(elem::Matrix{SymPy.Sym}, nodes::Array{Float64, 1};
    intNodes1 = 0, intNodes2 = 0, dof1 = 1, dof2 = 1)
```

## Assembling functions -- 2d

```@docs
SymFEL.assemble_squaremesh_FE_matrix(el_mat::Array{Float64, 2},
                                     elements::Array{Int64, 2};
                                     order1 = 1,
                                     order2 = 1,
                                     dof1 = 1,
                                     dof2 = 1)
```

```@docs
SymFEL.assemble_squaremesh_FE_matrix(el_mat::Array{Float64, 2},
                                     elements::Array{Int64, 2},
                                     el_labels::Array{Int64, 1};
                                     order1 = 1,
                                     order2 = 1,
                                     dof1 = 1,
                                     dof2 = 1)
```

```@docs
SymFEL.assemble1d_squaremesh_FE_matrix(el_mat::Array{Float64, 2},
                                       elements::Array{Int64, 2},
                                       elements1d::Array{Int64, 2};
                                       order1 = 1,
                                       order2 = 1,
                                       dof1 = 1,
                                       dof2 = 1)
```

## Assembling functions -- 3d

```@docs
SymFEL.assemble_cubemesh_FE_matrix(el_mat::Array{Float64, 2},
                                   elements::Array{UInt64, 2};
                                   order1 = 1,
                                   order2 = 1,
                                   dof1 = 1,
                                   dof2 = 1)
```

```@docs
SymFEL.assemble_cubemesh_FE_matrix(el_mat::Array{Float64, 2},
	                               elements::Array{UInt64, 2},
                                   el_labels::Array{UInt64, 1};
                                   order1 = 1,
                                   order2 = 1,
                                   dof1 = 1,
                                   dof2 = 1)
```

```@docs
SymFEL.assemble2d_cubemesh_FE_matrix(el_mat::Array{Float64, 2},
                                     elements::Array{UInt64, 2},
                                     elements2d::Array{UInt64, 2};
                                     order1 = 1,
                                     order2 = 1,
                                     dof1 = 1,
                                     dof2 = 1)
```
