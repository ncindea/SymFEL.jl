# SymFEL.jl Documentation


This package contains several functions usefull for the implementation of the finite elements method (FEM). We use the package `SymPy` for the computation of the basis finite element functions for both Lagrange and Hermite finite elements. We also compute the elementary matrices in 1d and for quadrangular elements in 2d. 
For the two-dimensional meshes we employ [Gmsh](https://gmsh.info/). The Gmsh SDK should be installed in your path following the instruction on the Gmsh site.

The utilization of the package is illustrated by several examples in 1d, 2d and 3d.

This package contains several functions usefull for the implementation of the finite elements method (FEM). We use the package `SymPy` for the computation of the basis finite element functions for both Lagrange and Hermite finite elements. We also compute the elementary matrices in 1d and for quadrangular elements in 2d and 3d obtained by tensor products from 1d elementary matrices. 
For the two and three dimensional meshes we employ [Gmsh](https://gmsh.info/). The Gmsh SDK should be installed in your path following the instruction on the Gmsh site.

The utilisation of the package is illustrated by several examples in 1d, 2d and 3d.


A documentation website is available at [http://lmbp.uca.fr/~cindea/software/SymFEL.jl/](http://lmbp.uca.fr/~cindea/software/SymFEL.jl/).

Here are some features to develop in the next versions:
- 2d triangular finite elements
- 3d tetrahedral finite elements
- combine the symbolic integration to numerical integration to speed up the execution.

## Contents
```@contents
Pages = ["index.md", "examples1.md", "examples2.md", "examples3.md"]
```

## SymFEl.jl

```@docs
SymFEL.get_em(deg1, deg2, der1::Function, der2::Function;
                fe1="Lagrange", fe2="Lagrange",
                x=symbols("x"), h=symbols("h"))
```

```@docs
SymFEL.get_em(deg1=1, deg2=1, der1=0, der2=0;
              fe1="Lagrange", fe2="Lagrange",
              x=symbols("x"), h=symbols("h"))
```

```@docs
SymFEL.get_etensor(deg1=1, deg2=1, deg3=1, der1=0, der2=0, der3=0;
                     fe1="Lagrange", fe2="Lagrange", fe3="Lagrange",
                     x=symbols("x"), h=symbols("h"))
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

```@docs
SymFEL.get_square_etensor(Tx::Array{SymPy.Sym, 3},
                          Ty::Array{SymPy.Sym, 3},
                          nc::Tuple{Array{Int64,1},Array{Int64,1}},
                          nr::Tuple{Array{Int64,1},Array{Int64,1}},
                          np::Tuple{Array{Int64,1},Array{Int64,1}})
```

## Lagrange finite elements

```@docs
SymFEL.get_lagrange_basis(n = 1, varcoeff = false;
	                      x = symbols("x"), h = symbols("h"))
```

```@docs
SymFEL.get_lagrange_em(p = 1, m = 0, n = 0;
                         x = symbols("x"), h = symbols("h"))
```

```@docs
SymFEL.get_square_lagrange_em((px, py) = (1, 1),
                                (mx, my) = (0, 0),
                                (nx, ny) = (0, 0);
                                x = symbols("x"),
                                hx = symbols("h"),
                                y = symbols("y"),
                                hy = symbols("h"))
```

```docs
SymFEL.get_lagrange_em_varcoeff(p = 1, m = 0, n = 0, f = 1;
                                  x = symbols("x"), h = symbols("h"),
                                  xa = symbols("xa"), xb = symbols("xb"))
```
## Hermite finite elements

```@docs
SymFEL.get_hermite_basis(n = 3, varcoeff = false;
                         x=symbols("x"), h=symbols("h"))
```

```@docs
SymFEL.get_hermite_em(p = 3, m = 0, n = 0;
                      x=symbols("x"), h=symbols("h"))
```

```@docs
SymFEL.get_hermite_em_varcoeff(p = 3, m = 0, n = 0, f = 1;
                               x=symbols("x"), h=symbols("h"))
```

```@doc
SymFEL.get_square_hermite_em((px, py) = (3, 3),
                               (mx, my) = (0, 0),
                               (nx, ny) = (0, 0);
                               x=symbols("x"), hx=symbols("hx"),
                               y=symbols("y"), hy=symbols("hy"))
```

```@docs
SymFEL.interpolate(fd, t, ti;
                   x=symbols("x"), h=symbols("h"))
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
SymFEL.assemble_1d_FE_matrix_multcoeff(elet::Array{Float64, 3},
                                coeff::Vector{Float64},
                                nbNodes::Int64;
                                intNodes1 = 0,
                                intNodes2 = 0,
                                intNodes3 = 0,
                                dof1 = 1,
                                dof2 = 1,
                                dof3 = 1) 
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
SymFEL.assemble_squaremesh_FE_matrix_coeffmult(el_ten::Array{Float64, 3},
                                                 coeff::Array{Float64, 1},
                                                 elements::Array{Int64, 2};
                                                 order1 = 1,
                                                 order2 = 1,
                                                 order3 = 1,
                                                 dof1 = 1,
                                                 dof2 = 1,
                                                 dof3 = 1)
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

## Exports elementary matrices and tensors to pure julia code

```@docs
SymFEL.exports_mat(fun::String, mat::Matrix{Sym}; path::String="./")
```

```@docs
SymFEL.exports_ten(fun::String, ten::Array{Sym, 3}; path="./")
```
