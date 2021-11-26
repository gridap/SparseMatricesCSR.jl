module SparseMatricesCSR

using SparseArrays
using LinearAlgebra
using SuiteSparse

import Base: convert, copy, size, getindex, setindex!, show, count, *, IndexStyle
import LinearAlgebra: mul!, lu, lu!
import SparseArrays: nnz, getnzval, nonzeros, nzrange
import SparseArrays: findnz, rowvals, getnzval, issparse

export SparseMatrixCSR
export SymSparseMatrixCSR
export sparsecsr, symsparsecsr
export colvals, getBi, getoffset

include("SparseMatrixCSR.jl")

include("SymSparseMatrixCSR.jl")

end # module

# DONE sparse
# DONE issparse
# DONE nnz
# DONE nonzeros
# DONE colvals
# DONE nzrange
# DONE findnz
# DONE getrowptr
# DONE getnzval
# DONE getcolval
# TODO spzeros
# TODO spdiagm
# TODO blockdiag
# TODO sprand
# TODO sprandn
# TODO droptol!
# TODO dropzeros
# TODO dropzeros!
# TODO permute
# TODO permute!
