#!/bin/sh
julia -e 'include("docs/make.jl")'
scp -P 4220 -r docs/build/* localhost:~/.public_html/software/FEMTools.jl
