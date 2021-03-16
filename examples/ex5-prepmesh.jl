import gmsh
using JLD2

gmsh.initialize()
gmsh.open("cube.msh")
nodes = gmsh.model.mesh.getNodes()
nodes_label = nodes[1]
nodes_N = length(nodes_label)
nodes_coordinate = reshape(nodes[2], (3, nodes_N))
nodes_boundary = unique(gmsh.model.mesh.getNodesByElementType(3)[1])
nodes_boundary_N = length(nodes_boundary)

elements = gmsh.model.mesh.getElements()

elements_bound_label = elements[2][1]
elements_bound_N = length(elements_bound_label)

elements_int_label = elements[2][2]
elements_int_N = length(elements_int_label)

elements_bound = reshape(elements[3][1], (4, elements_bound_N))
elements_int = reshape(elements[3][2], (8, elements_int_N))
gmsh.finalize()
@save "cube.jld2"
