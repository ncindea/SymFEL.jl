# todo : document and odd some other nasty functions
"""
    assemble_1d_FE_matrix(elem, nbNodes)
"""
function assemble_1d_FE_matrix(elem, nbNodes)
  M = spzeros(Float64, nbNodes, nbNodes)
  for i = 1:nbNodes - 1
    M[i:(i + 1), i:(i + 1)] = M[i:(i + 1), i:(i + 1)] + elem
  end
  M
end
