module SparseMatricesCSR

using SparseArrays
using LinearAlgebra

import Base: size, getindex, show, count
import SparseArrays: nnz, getnzval, nzvalview, nonzeros, nzrange, findnz

export SparseMatrixCSR
export SymSparseMatrixCSR
export push_coo!, finalize_coo!, sparsecsr, symsparsecsr, getrowptr, getcolval

include("SparseMatrixCSC.jl")
include("SparseMatrixCSR.jl")
include("SymSparseMatrixCSR.jl")

end # module
