using SparseArrays

"""
    assemble_cubemesh_FE_matrix_coeffmult(el_ten::Array{Float64, 3},
                                  coeff::Array{Float64, 1}, 
                                  elements::Array{Int64, 2};
                                  order1 = 1, order2 = 1, order3 = 1
                                  dof1 = 1, dof2 = 1, dof3 = 1)

Assemble a finite elements matrix corresponding to a 2 dimensional square mesh.

# Arguments
  * `el_ten`   : elementary finite elements matrix
  * `coeff`    : variable coefficient
  * `elements` : list of elements
  * `order1`   : order for lhs
  * `order2`   : order for rhs
  * `order3`   : order for phs
  * `dof1`     : number of degrees of freedom for each node for lhs
  * `dof2`     : number of degrees of freedom for each node for rhs
  * `dof3`     : number of degrees of freedom for each node for phs
"""
function assemble_cubemesh_FE_matrix_coeffmult(el_ten::Array{Float64, 3},
                                               coeff::Array{Float64, 1},
                                               elements::Array{Int64, 2};
                                               order1 = 1,
                                               order2 = 1,
                                               order3 = 1,
                                               dof1 = 1,
                                               dof2 = 1,
                                               dof3 = 1)
    
    n_order1 = (order1 + 1)^3
    n_order2 = (order2 + 1)^3
    n_order3 = (order3 + 1)^3

    el1 = copy(elements)
    el2 = copy(elements)
    el3 = copy(elements)

    if n_order1 == 1
        nodes1 = 1:size(elements, 2)
        el1[1,:] = nodes1
    else
        nodes1 = sort(unique(elements[1:n_order1,:][:]))
    end
    if n_order2 == 1
        nodes2 = 1:size(elements, 2)
        el2[1,:] = nodes2
    else
        nodes2 = sort(unique(elements[1:n_order2,:][:]))
    end
    if n_order3 == 1
        nodes3 = 1:size(elements, 2)
        el3[1,:] = nodes3
    else
        nodes3 = sort(unique(elements[1:n_order3,:][:]))
    end
   
    nodes1_N = length(nodes1)
    nodes2_N = length(nodes2)
    nodes3_N = length(nodes3)

    elements_N = size(elements, 2)

    M = spzeros(Float64, nodes1_N * dof1, nodes2_N * dof2)

    ls = size(el_ten, 1)
    rs = size(el_ten, 2)
    ps = size(el_ten, 3)
    l1 = zeros(UInt64, n_order1 * dof1)
    l2 = zeros(UInt64, n_order2 * dof2)
    l3 = zeros(UInt64, n_order3 * dof3)
    v1 = @. dof1 * ((1:n_order1) - 1)
    v2 = @. dof2 * ((1:n_order2) - 1)
    v3 = @. dof3 * ((1:n_order3) - 1)

    for i = 1:elements_N
        el_mat = zeros(ls, rs)
        for j = 1:dof1
            l1[v1 .+ j] = @. dof1 * (el1[1:n_order1, i] - 1) + j
        end
        for j = 1:dof2
            l2[v2 .+ j] = @. dof2 * (el2[1:n_order2, i] - 1) + j
        end
        for j = 1:dof3
            l3[v3 .+ j] = @. dof3 * (el3[1:n_order3, i] - 1) + j
        end
        for ip = 1:ps
            el_mat += coeff[l3[ip]] * el_ten[:, :, ip]
        end
        
        M[l1, l2] +=  el_mat
    end
    M
end


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
    n_order1 = (order1+1)^3
    n_order2 = (order2+1)^3

    el1 = copy(elements)
    el2 = copy(elements)

    if n_order1 == 1
        nodes1 = 1:size(elements, 2)
        el1[1,:] = nodes1
    else
        nodes1 = sort(unique(elements[1:n_order1,:][:]))
    end
    if n_order2 == 1
        nodes2 = 1:size(elements, 2)
        el2[1,:] = nodes2
    else
        nodes2 = sort(unique(elements[1:n_order2,:][:]))
    end
    
   
    nodes1_N = length(nodes1)
    nodes2_N = length(nodes2)

    elements_N = size(elements, 2)
    M = spzeros(Float64, nodes1_N * dof1, nodes2_N * dof2)
    
    l1 = zeros(Int64, n_order1 * dof1)
    l2 = zeros(Int64, n_order2 * dof2)

    v1 = @. dof1 * ((1:n_order1)-1)
    v2 = @. dof2 * ((1:n_order2)-1)
    
    for i = 1:elements_N
        for j = 1:dof1
            l1[v1  .+ j] = @. dof1 * (el1[1:n_order1, i] - 1) + j
        end
        for j = 1:dof2
            l2[v2  .+ j] = @. dof2 * (el2[1:n_order2, i] - 1) + j
        end
        
        M[l1, l2] += el_mat
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
    
    n_order1 = (order1+1)^3
    n_order2 = (order2+1)^3
    
    el1 = copy(elements)
    el2 = copy(elements)

    if n_order1 == 1
        nodes1 = 1:size(elements, 2)
        el1[1,:] = nodes1
    else
        nodes1 = sort(unique(elements[1:n_order1,:][:]))
    end
    if n_order2 == 1
        nodes2 = 1:size(elements, 2)
        el2[1,:] = nodes2
    else
        nodes2 = sort(unique(elements[1:n_order2,:][:]))
    end
    
    
    nodes1_N = length(nodes1)
    nodes2_N = length(nodes2)

    elements_N = size(elements, 2)
    labels_N = length(el_labels)
    M = spzeros(Float64, nodes1_N * dof1, nodes2_N * dof2)
    
    l2 = zeros(Int64, n_order1 * dof1)
    r2 = zeros(Int64, n_order2 * dof2)

    v1 = @. dof1 * ((1:n_order1)-1)
    v2 = @. dof2 * ((1:n_order2)-1)

    
    for i = 1:labels_N
        for j = 1:dof1
            l2[v1 .+ j] = @. dof1 * (el1[1:n_order1, el_labels[i]] - 1) + j
        end
        for j = 1:dof2
            r2[v2 .+ j] = @. dof2 * (el2[1:n_order2, el_labels[i]] - 1) + j
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

    n_order1 = (order1+1)^3
    n_order2 = (order2+1)^2

    nodes = sort(unique(elements[1:n_order1,:][:]))

    nodes_N = length(nodes)
    

    elements2d_N = size(elements2d, 2)
    M = spzeros(Float64, nodes_N * dof1, nodes_N * dof2)
    l2 = zeros(Int64, n_order2 * dof2)
    r2 = zeros(Int64, n_order2 * dof2)
    v1 = dof2 * ((1:n_order2).-1)
    v2 = dof2 * ((1:n_order2).-1)
    
    for i = 1:elements2d_N
        for j = 1:dof2
            l2[v1 .+ j] = @. dof1 * (elements2d[1:n_order2, i] - 1) + j
        end
        for j = 1:dof2
            r2[v2 .+ j] = @. dof2 * (elements2d[1:n_order2, i] - 1) + j
        end
        M[l2, r2] +=  el_mat
    end
    M
end
