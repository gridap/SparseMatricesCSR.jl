module SparseMatricesCSR

using SparseArrays
using LinearAlgebra

import Base: size, getindex, show, count
import SparseArrays: nnz, getnzval, nzvalview, nonzeros, nzrange, findnz

export SparseMatrixCSR
export SymSparseMatrixCSR
export sparsecsr, symsparsecsr, getrowptr, getcolval

include("SparseMatrixCSR.jl")
include("SymSparseMatrixCSR.jl")

end # module
