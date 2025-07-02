# todo : document and add some other nasty functions
using SparseArrays

"""
    Mesh1d

A type describing a mesh for a one dimensional domain.
"""
struct Mesh1d
    nodes::Array{Float64, 1}
    elements::Array{Int64, 2}
end

"""
    function assemble_1d_FE_matrix(elem::Array{Float64, 2}, nbNodes::Int64;
                                   intNodes1 = 0, intNodes2 = 0, dof1 = 1, dof2 = 1)

Assemble a finite elements matrix corresponding to a 1 dimensional uniform mesh.

# Arguments
  * `elem` : elementary finite elements matrix
  * `nbNodes` : number of nodes in the mesh
  * `intNodes1` : number of interior nodes in the interior of each element for lhs
  * `intNodes2` : number of interior nodes in the interior of each element for rhs
  * `dof1`     : number of degrees of freedom for each node for lhs
  * `dof2`     : number of degrees of freedom for each node for rhs
"""
function assemble_1d_FE_matrix(elem::Array{Float64, 2}, nbNodes::Int64;
                               intNodes1 = 0, intNodes2 = 0, dof1 = 1, dof2 = 1)

    nbNodesTotal1 = (nbNodes + (nbNodes - 1) * intNodes1)
    nbNodesTotal2 = (nbNodes + (nbNodes - 1) * intNodes2)

    M = spzeros(Float64, nbNodesTotal1 * dof1, nbNodesTotal2 * dof2)
    for i = 1:(nbNodes - 1)
        l1 = (i - 1 + (i - 1) * intNodes1) * dof1 + 1
        r1 = (i + 1 + i * intNodes1) * dof1
        l2 = (i - 1 + (i - 1) * intNodes2) * dof2 + 1
        r2 = (i + 1 + i * intNodes2) * dof2
        M[l1:r1, l2:r2] += elem
    end

    return M
end


"""
    assemble_1d_FE_matrix_multcoeff(elet::Array{Float64, 3},
                                    coeff::Vector{Float64},
                                    nbNodes::Int64;
                                    intNodes1 = 0,
                                    intNodes2 = 0,
                                    intNodes3 = 0,
                                    dof1 = 1,
                                    dof2 = 1,
                                    dof3 = 1)
Assemble a finite elements matrix corresponding to a 1 dimensional uniform mesh with a variable coefficient.

# Arguments
  * `elet` : elementary finite elements tensor
  * `coeff` : multiplicative coefficient (corresponding to 3rd FEM)
  * `nbNodes` : number of nodes in the mesh
  * `intNodes1` : number of interior nodes in the interior of each element for lhs
  * `intNodes2` : number of interior nodes in the interior of each element for rhs
  * `intNodes4` : number of interior nodes in the interior of each element for phs
  * `dof1`     : number of degrees of freedom for each node for lhs
  * `dof2`     : number of degrees of freedom for each node for rhs
  * `dof3`     : number of degrees of freedom for each node for phs
"""
function assemble_1d_FE_matrix_multcoeff(elet::Array{Float64, 3},
                                         coeff::Vector{Float64},
                                         nbNodes::Int64;
                                         intNodes1 = 0,
                                         intNodes2 = 0,
                                         intNodes3 = 0,
                                         dof1 = 1,
                                         dof2 = 1,
                                         dof3 = 1)

    nbNodesTotal1 = (nbNodes + (nbNodes - 1) * intNodes1)
    nbNodesTotal2 = (nbNodes + (nbNodes - 1) * intNodes2)
    ls = size(elet, 1)
    rs = size(elet, 2)
    ps = size(elet, 3)


    M = spzeros(Float64, nbNodesTotal1 * dof1, nbNodesTotal2 * dof2)

    for i = 1:(nbNodes - 1)
        elem = zeros(Float64, ls, rs)
        
        l1 = (i - 1 + (i - 1) * intNodes1) * dof1 + 1
        r1 = (i + 1 + i * intNodes1) * dof1
        l2 = (i - 1 + (i - 1) * intNodes2) * dof2 + 1
        r2 = (i + 1 + i * intNodes2) * dof2
        l3 = (i - 1 + (i - 1) * intNodes3) * dof3 + 1
        r3 = (i + 1 + i * intNodes3) * dof3
        for ip = l3:r3
            elem += coeff[ip] * elet[:,:,ip - l3 + 1]
        end
        
        M[l1:r2, l2:r2] += elem
    end
    M
end

"""
    assemble_1d_FE_symmatrix(elem::Array{Sym, 2}, nbNodes::Int64;
      intNodes = 0, dof1 = 1, dof2 = 1)

Assemble a symbolic finite elements matrix corresponding to a 1 dimensional uniform mesh.

# Arguments
  * `elem` : elementary finite elements matrix
  * `nbNodes` : number of nodes in the mesh
  * `intNodes1` : number of interior nodes in the interior of each element for lhs
  * `intNodes2` : number of interior nodes in the interior of each element for rhs
  * `dof1`     : number of degrees of freedom for each node for lhs
  * `dof2`     : number of degrees of freedom for each node for rhs
"""
function assemble_1d_FE_symmatrix(elem::Array{Sym, 2}, nbNodes::Int64;
                                  intNodes1 = 0, intNodes2 = 0, dof1 = 1, dof2 = 1)

    nbNodesTotal1 = (nbNodes + (nbNodes - 1) * intNodes1)
    nbNodesTotal2 = (nbNodes + (nbNodes - 1) * intNodes2)

    M = spzeros(Float64, nbNodesTotal1 * dof1, nbNodesTotal2 * dof2)
    for i = 1:(nbNodes - 1)
        l1 = (i - 1 + (i - 1) * intNodes1) * dof1 + 1
        r1 = (i + 1 + i * intNodes1) * dof1
        l2 = (i - 1 + (i - 1) * intNodes2) * dof2 + 1
        r2 = (i + 1 + i * intNodes2) * dof2
        M[l1:r1, l2:r2] += elem
    end
    M
end


"""
    assemble_1d_nu_FE_matrix_varcoeff(elem::Matrix{SymPy.Sym}, nodes::Array{Float64, 1};
      intNodes = 0, dof1 = 1, dof2 = 1)

Assemble a finite elements matrix corresponding to a 1 dimensional non-uniform mesh.

# Arguments
  * `elem` : elementary finite elements matrix (with elements of type SymPy.Sym)
  * `nodes` : vector of nodes composing the non uniform-mesh
  * `intNodes1` : number of interior nodes in the interior of each element for lhs
  * `intNodes2` : number of interior nodes in the interior of each element for rhs
  * `dof1`     : number of degrees of freedom for each node for lhs
  * `dof2`     : number of degrees of freedom for each node for rhs
"""
function assemble_1d_nu_FE_matrix_varcoeff(elem::Matrix{SymPy.Sym},
                                           nodes::Array{Float64, 1};
                                           intNodes1 = 0,
                                           intNodes2 = 0, dof1 = 1, dof2 = 1,
                                           xa=symbols("xa"), xb=symbols("xb"))

    nbNodes = length(nodes)
    nbNodesTotal1 = (nbNodes + (nbNodes - 1) * intNodes1)
    nbNodesTotal2 = (nbNodes + (nbNodes - 1) * intNodes2)

    elem_loc = lambdify(elem, (xa, xb), use_julia_code=true)

    M = spzeros(Float64, nbNodesTotal1 * dof1, nbNodesTotal2 * dof2)
    for i = 1:(nbNodes - 1)
        l1 = (i - 1 + (i - 1) * intNodes1) * dof1 + 1
        r1 = (i + 1 + i * intNodes1) * dof1
        l2 = (i - 1 + (i - 1) * intNodes2) * dof2 + 1
        r2 = (i + 1 + i * intNodes2) * dof2
        M[l1:r1, l2:r2] += convert(Matrix{Float64}, elem_loc(nodes[i], nodes[i+1]))
    end 
    M
end

"""
    assemble_1d_nu_FE_matrix(elem::Matrix{SymPy.Sym}, nodes::Array{Float64, 1};
      intNodes = 0, dof1 = 1, dof2 = 1)

Assemble a finite elements matrix corresponding to a 1 dimensional non-uniform mesh.

# Arguments
  * `elem` : elementary finite elements matrix (with elements of type SymPy.Sym)
  * `nodes` : vector of nodes composing the non uniform-mesh
  * `intNodes1` : number of interior nodes in the interior of each element for lhs
  * `intNodes2` : number of interior nodes in the interior of each element for rhs
  * `dof1`     : number of degrees of freedom for each node for lhs
  * `dof2`     : number of degrees of freedom for each node for rhs
"""
function assemble_1d_nu_FE_matrix(elem::Matrix{SymPy.Sym}, nodes::Array{Float64, 1};
                                  intNodes1 = 0, intNodes2 = 0, dof1 = 1, dof2 = 1, h=symbols("h"))

    nbNodes = length(nodes)
    nbNodesTotal1 = (nbNodes + (nbNodes - 1) * intNodes1)
    nbNodesTotal2 = (nbNodes + (nbNodes - 1) * intNodes2)

    M = spzeros(Float64, nbNodesTotal1 * dof1, nbNodesTotal2 * dof2)
    elem_loc = lambdify(elem, use_julia_code=true)
    
    for i = 1:(nbNodes - 1)
        l1 = (i - 1 + (i - 1) * intNodes1) * dof1 + 1
        r1 = (i + 1 + i * intNodes1) * dof1
        l2 = (i - 1 + (i - 1) * intNodes2) * dof2 + 1
        r2 = (i + 1 + i * intNodes2) * dof2
        M[l1:r1, l2:r2] += convert(Matrix{Float64}, elem_loc(nodes[i + 1] - nodes[i]))
    end
    M
end


"""
   as1d_FEM is an alias for assemble_1d_FE_matrix
"""
as1d_FEM = assemble_1d_FE_matrix

"""
   as1d_nuFEM is an alias for assemble_1d_nu_FE_matrix
"""
as1d_nuFEM = assemble_1d_nu_FE_matrix

"""
   as1d_nuFEM is an alias for assemble_1d_FE_matrix_multcoeff
"""
as1d_FEM_mc = assemble_1d_FE_matrix_multcoeff
