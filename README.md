# SymFEL.jl


| Testing                                                                                                       |
|---------------------------------------------------------------------------------------------------------------|
| [![Build Status](https://travis-ci.org/ncindea/SymFEL.jl.svg?branch=master)](https://travis-ci.org/ncindea/SymFEL.jl) |


This package contains several functions usefull for the implementation of the finite elements method (FEM). We use the package `SymPy` for the computation of the basis finite element functions for both Lagrange and Hermite finite elements. We also compute the elementary matrices in 1d and for quadragular elements in 2d and 3d. 
For the two and three dimensional meshes we employ [Gmsh](https://gmsh.info/). The Gmsh SDK should be installed in your path following the instruction on the Gmsh site.

The utilisation of the package is illustrated by several examples in 1d, 2d and 3d.

A documentation website is available at [http://lmbp.uca.fr/~cindea/software/SymFEL.jl/](http://lmbp.uca.fr/~cindea/software/SymFEL.jl/).

Here are some features to develop in the next versions:
- 2d triangular finite elements
- 3d tetrahedral finite elements
- combine the symbolic integration to numerical integration to speed-up the execution
- add more high-level functionalities.
