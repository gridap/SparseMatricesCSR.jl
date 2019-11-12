module SparseMatricesCSR

using SparseArrays
using LinearAlgebra

import Base: size, getindex, show, count
import LinearAlgebra: mul!
import SparseArrays: nnz, getnzval, nzvalview, nonzeros, nzrange, findnz, rowvals

export SparseMatrixCSR
export SymSparseMatrixCSR
export push_coo!, finalize_coo!, sparsecsr, symsparsecsr, colvals
export hasrowmajororder, hascolmajororder, getptr, getindices

include("SparseMatrixCSC.jl")
include("SparseMatrixCSR.jl")
include("SymSparseMatrixCSR.jl")

end # module
