module SparseMatrix

using SparseArrays

import Base: size, getindex, show, count
import SparseArrays: nnz, getnzval, nzvalview, nonzeros, nzrange, findnz

export SparseMatrixCSR
export sparsecsr, getrowptr, getcolval

include("SparseMatrixCSR.jl")

end # module
