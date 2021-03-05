using SparseArrays

"""
    assemble_cubemesh_FE_matrix(el_mat::Array{Float64, 2},
                                elements::Array{UInt64, 2};
                                order1 = 1, order2 = 1,
                                dof1 = 1, dof2 = 1)

Assemble a finite elements matrix corresponding to a 3 dimensional cube mesh.

# Arguments
  * `el_mat`   : elementary finite elements matrix
  * `elements` : list of elements
  * `order1`   : order for lhs
  * `order2`   : order for rhs
  * `dof1`     : number of degrees of freedom for each node for lhs
  * `dof2`     : number of degrees of freedom for each node for rhs
"""
function assemble_cubemesh_FE_matrix(el_mat::Array{Float64, 2},
                                     elements::Array{UInt64, 2};
                                     order1 = 1,
                                     order2 = 1,
                                     dof1 = 1,
                                     dof2 = 1)
    
    n_order1 = (order1+1)*(4 + (order1 - 1) * (order1 + 3))
    n_order2 = (order2+1)*(4 + (order2 - 1) * (order2 + 3))
    nodes1 = sort(unique(elements[1:n_order1,:][:]))
    nodes2 = sort(unique(elements[1:n_order2,:][:]))
    nodes1_N = length(nodes1)
    nodes2_N = length(nodes2)
    
    nodes1_i = zeros(Int64, maximum(nodes1))
    nodes1_i[nodes1] = 1:nodes1_N
    nodes2_i = zeros(Int64, maximum(nodes2))
    nodes2_i[nodes2] = 1:nodes2_N

    elements_N = size(elements, 2)
    M = spzeros(Float64, nodes1_N * dof1, nodes2_N * dof2)
    l2 = zeros(Int64, n_order1 * dof1)
    r2 = zeros(Int64, n_order2 * dof2)
    for i = 1:elements_N
        l = elements[1:n_order1, i]
        r = elements[1:n_order2, i]
        l = nodes1_i[l]
        r = nodes2_i[r]
        for j = 1:dof1
            l2[dof1 * ((1:n_order1).-1) .+ j] = dof1 * (l .- 1)  .+ j
        end
        for j = 1:dof2
            r2[dof2 * ((1:n_order2).-1) .+ j] = dof2 * (r .- 1)  .+ j
        end
        M[l2, r2] +=  el_mat
    end
    M
end


"""
    assemble_cubemesh_FE_matrix(el_mat::Array{Float64, 2},
                                  elements::Array{UInt64, 2},
                                  el_labels::Array{UInt64, 1};
                                  order1 = 1, 
                                  order2 = 1,
                                  dof1 = 1, 
                                  dof2 = 1)

Assemble a finite elements matrix corresponding to a 3 dimensional cube mesh.

# Arguments
  * `el_mat`   : elementary finite elements matrix
  * `elements` : list of elements
  * `el_labels`: list of element numbers on which we integrate
  * `order1`   : order for lhs
  * `order2`   : order for rhs
  * `dof1`     : number of degrees of freedom for each node for lhs
  * `dof2`     : number of degrees of freedom for each node for rhs
"""
function assemble_cubemesh_FE_matrix(el_mat::Array{Float64, 2},
                                     elements::Array{UInt64, 2},
                                     el_labels::Array{UInt64, 1};
                                     order1 = 1,
                                     order2 = 1,
                                     dof1 = 1,
                                     dof2 = 1)
    
    n_order1 = (order1+1)*(4 + (order1 - 1) * (order1 + 3))
    n_order2 = (order2+1)*(4 + (order2 - 1) * (order2 + 3))
    nodes1 = sort(unique(elements[1:n_order1,:][:]))
    nodes2 = sort(unique(elements[1:n_order2,:][:]))
    nodes1_N = length(nodes1)
    nodes2_N = length(nodes2)
    
    nodes1_i = zeros(Int64, maximum(nodes1))
    nodes1_i[nodes1] = 1:nodes1_N
    nodes2_i = zeros(Int64, maximum(nodes2))
    nodes2_i[nodes2] = 1:nodes2_N

    elements_N = size(elements, 2)
    labels_N = length(el_labels)
    M = spzeros(Float64, nodes1_N * dof1, nodes2_N * dof2)
    l2 = zeros(Int64, n_order1 * dof1)
    r2 = zeros(Int64, n_order2 * dof2)
    for i = 1:labels_N
        l = elements[1:n_order1, el_labels[i]]
        r = elements[1:n_order2, el_labels[i]]
        l = nodes1_i[l]
        r = nodes2_i[r]
        for j = 1:dof1
            l2[dof1 * ((1:n_order1).-1) .+ j] = dof1 * (l .- 1)  .+ j
        end
        for j = 1:dof2
            r2[dof2 * ((1:n_order2).-1) .+ j] = dof2 * (r .- 1)  .+ j
        end
        M[l2, r2] +=  el_mat
    end
    M
end

"""
    assemble2d_cubemesh_FE_matrix(el_mat::Array{Float64, 2},
                                  elements::Array{UInt64, 2},
                                  el_labels::Array{UInt64, 1};
                                  order1 = 1, 
                                  order2 = 1,
                                  dof1 = 1, 
                                  dof2 = 1)

Assemble a finite elements matrix corresponding to a 3 dimensional cube mesh.

# Arguments
  * `el_mat`   : elementary finite elements matrix
  * `elements` : list of elements
  * `elements1s`: list of 1d elements
  * `order1`   : order for lhs
  * `order2`   : order for rhs
  * `dof1`     : number of degrees of freedom for each node for lhs
  * `dof2`     : number of degrees of freedom for each node for rhs
"""
function assemble2d_cubemesh_FE_matrix(el_mat::Array{Float64, 2},
                                       elements::Array{UInt64, 2},
                                       elements2d::Array{UInt64, 2};
                                       order1 = 1,
                                       order2 = 1,
                                       dof1 = 1,
                                       dof2 = 1)
    
    n_order1 = (order1+1)*(4 + (order1 - 1) * (order1 + 3))
    n_order2 = 4 + (order1 - 1) * (order1 + 3)

    nodes = sort(unique(elements[1:n_order1,:][:]))

    nodes_N = length(nodes)
    

    elements2d_N = size(elements2d, 2)
    M = spzeros(Float64, nodes_N * dof1, nodes_N * dof2)

    r2 = zeros(Int64, n_order2 * dof2)
    l2 = zeros(Int64, n_order2 * dof2)
    for i = 1:elements2d_N
        l = elements1d[1:n_order2, i]
        r = elements1d[1:n_order2, i]
        
        
        for j = 1:dof2
            l2[dof2 * ((1:n_order2).-1) .+ j] = dof1 * (l .- 1)  .+ j
        end
        for j = 1:dof2
            r2[dof2 * ((1:n_order2).-1) .+ j] = dof2 * (r .- 1)  .+ j
        end
        M[l2, r2] +=  el_mat
    end
    M
end
