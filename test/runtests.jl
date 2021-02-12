using SparseArrays
using LinearAlgebra
using SparseMatricesCSR
using Test

@testset "SparseMatrixCSR" begin include("SparseMatrixCSR.jl") end

@testset "SymSparseMatrixCSR" begin include("SymSparseMatrixCSR.jl") end
