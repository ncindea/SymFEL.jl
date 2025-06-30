using SymPy
using JuliaFormatter

"""
   function exports_mat(fun::String, mat::Matrix{Sym}; path="./")

Exports a matrix of type `Matrix{Sym}` to a function written in julia.
The function named `fun` is saved in the file `strin(path, fun, ".jl")`.
"""
function exports_mat(fun::String, mat::Matrix{Sym}; path::String="./")
    file_name = string(expanduser(path), fun, ".jl")
    
    s = Set([])
    for i = 1:size(mat, 1)
        for j = 1:size(mat, 2)
            union!(s, free_symbols(mat[i, j]))
        end
    end
    
    nums = length(s)
    svec = String[]
    
    for e in s
        push!(svec, SymPy.sympy.julia_code(e))
    end
    sort!(svec)
    head = string("function ", fun, "(")
    for n = 1:nums-1
        head = string(head, svec[n], "::Float64, ")
    end
    head = string(head, svec[nums], "::Float64)")

    body = SymPy.sympy.julia_code(mat)
    io = open(file_name, "w");
    write(io, "# This function was generated with 💙 by SymFEL.exports_mat().\n")
    write(io, head)
    write(io, "\n")
    write(io, body)
    write(io, "\n")
    write(io, "end")
    close(io)
    format_file(file_name)
end

"""
   function exports_ten(fun::String, ten::Array{Sym, 3}; path="./")

Exports a tensor of type `Array{Sym, 3}` to a function written in julia.
The function named `fun` is saved in the file `strin(path, fun, ".jl")`.
"""
function exports_ten(fun::String, ten::Array{Sym, 3}; path="./")
    file_name = string(expanduser(path), fun, ".jl")
    
    s = Set([])
    for i = 1:size(ten, 1)
        for j = 1:size(ten, 2)
            for k = 1:size(ten, 3)
                union!(s, free_symbols(ten[i, j, k]))
            end
            
        end
    end
    
    nums = length(s)
    svec = String[]
    
    for e in s
        push!(svec, SymPy.sympy.julia_code(e))
    end
    sort!(svec)
    head = string("function ", fun, "(")
    for n = 1:nums-1
        head = string(head, svec[n], "::Float64, ")
    end
    head = string(head, svec[nums], "::Float64)")

    body = String[]
    for n = 1:size(ten, 3)
        push!(body, SymPy.sympy.julia_code(ten[:, :, n]))
    end
    
    io = open(file_name, "w");
    write(io, "# This function was generated with 💙 by SymFEL.exports_ten().\n")
    write(io, head)
    write(io, "\n")
    write(io, "ten = zeros(Float64, ", string(size(ten)), ")\n")
    for n = 1:size(ten, 3)
        write(io, "ten[:, :, ", string(n), "] = ", body[n], "\n")
    end
    
    write(io, "\n")
    write(io, "ten\n")
    write(io, "end")
    close(io)
    format_file(file_name)
end
