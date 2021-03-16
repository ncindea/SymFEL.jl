import gmsh
using JLD2

gmsh.initialize()
gmsh.open("square-simple.msh")
nodes = gmsh.model.mesh.getNodes()
nodes_label = convert(Array{Int64}, nodes[1])
nodes_N = length(nodes_label)
nodes_coordinate = reshape(nodes[2], (3, nodes_N))
nodes_boundary = convert(Array{Int64}, unique(gmsh.model.mesh.getNodesByElementType(1)[1]))
nodes_boundary_N = length(nodes_boundary)

elements = gmsh.model.mesh.getElements()

elements_bound_label = convert(Array{Int64}, elements[2][1])
elements_bound_N = length(elements_bound_label)

elements_int_label = convert(Array{Int64}, elements[2][2])
elements_int_N = length(elements_int_label)

elements_bound = reshape(convert(Array{Int64}, elements[3][1]), (2, elements_bound_N))
elements_int = reshape(convert(Array{Int64}, elements[3][2]), (4, elements_int_N))
#gmsh.fltk.run()
gmsh.finalize()
@save "square-simple.jld2"
