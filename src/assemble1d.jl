# todo : document and add some other nasty functions
using SparseArrays
"""
    Mesh1d

A type describing a mesh for a one dimensional domain
"""
struct Mesh1d
  nodes::Array{Float64, 1}
  elements::Array{Int64, 2}
end

"""
    assemble_1d_FE_matrix(elem::Array{Float64, 2}, nbNodes::Int64;
      intNodes = 0, dof1 = 1, dof2 = 1)

Assemble a finite elements matrix corresponding to a 1 dimensional uniform mesh.

# Arguments
  * `elem` : elementary finite elements matrix
  * `nbNodes` : number of nodes in the mesh
  * `intNodes` : number of interior nodes in the interior of each elements
  * `dof1`     : number of degrees of freedom for each node for lhs
  * `dof2`     : number of degrees of freedom for each node for rhs
"""
function assemble_1d_FE_matrix(elem::Array{Float64, 2}, nbNodes::Int64;
    intNodes = 0, dof1 = 1, dof2 = 1)

  nbNodesTotal = (nbNodes + (nbNodes - 1) * intNodes)

  M = spzeros(Float64, nbNodesTotal * dof1, nbNodesTotal * dof2)
  for i = 1:nbNodes - 1
    l1 = (i - 1 + (i - 1) * intNodes) * dof1 + 1
    r1 = (i + 1 + i * intNodes) * dof1
    l2 = (i - 1 + (i - 1) * intNodes) * dof2 + 1
    r2 = (i + 1 + i * intNodes) * dof2
    M[l1:r1, l2:r2] = M[l1:r1, l2:r2] + elem
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
  * `intNodes` : number of interior nodes in the interior of each elements
  * `dof1`     : number of degrees of freedom for each node for lhs
  * `dof2`     : number of degrees of freedom for each node for rhs
"""
function assemble_1d_nu_FE_matrix(elem::Matrix{SymPy.Sym}, nodes::Array{Float64, 1};
    intNodes = 0, dof1 = 1, dof2 = 1)
  h = SymPy.symbols("h")
  nbNodes = length(nodes)
  nbNodesTotal = (nbNodes + (nbNodes - 1) * intNodes)

  M = spzeros(Float64, nbNodesTotal * dof1, nbNodesTotal * dof2)
  for i = 1:nbNodes - 1
    l1 = (i - 1 + (i - 1) * intNodes) * dof1 + 1
    r1 = (i + 1 + i * intNodes) * dof1
    l2 = (i - 1 + (i - 1) * intNodes) * dof2 + 1
    r2 = (i + 1 + i * intNodes) * dof2
    elem_loc = elem.subs(h, nodes[i + 1] - nodes[i])
    M[l1:r1, l2:r2] = M[l1:r1, l2:r2] + elem_loc
  end
  M
end
