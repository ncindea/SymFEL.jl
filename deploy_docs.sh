#!/bin/sh
julia --color=yes -e 'include("docs/make.jl")'
rsync -rvz -e 'ssh -p 4220' --progress docs/build/ cindea@localhost:/u2/users/projets/cindea/.public_html/software/SymFEL.jl

