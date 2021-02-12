using SparseArrays
using LinearAlgebra
using SparseMatricesCSR
using Test

#include("SparseMatrixCSC.jl")
@testset "SparseMatrixCSR" begin include("SparseMatrixCSR.jl") end
#include("SymSparseMatrixCSR.jl")
