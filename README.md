# SparseMatricesCSR

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gridap.github.io/SparseMatricesCSR.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gridap.github.io/SparseMatricesCSR.jl/dev)
[![Build Status](https://travis-ci.com/gridap/SparseMatricesCSR.jl.svg?branch=master)](https://travis-ci.com/gridap/SparseMatricesCSR.jl)
[![Codecov](https://codecov.io/gh/gridap/SparseMatricesCSR.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gridap/SparseMatricesCSR.jl)

Compressed Sparse Row (CSR) sparse matrix implementations.

It includes `SparseMatrixCSR` and `SymSparseMatrixCSR` wrappers.

## Basic Usage

```julia
julia> push_coo!(SymSparseMatrixCSR,I,J,V,i1,j1,v1)
julia> ...
julia> push_coo!(SymSparseMatrixCSR,I,J,V,in,jn,vn)
julia> finalize_coo!(SparseMatrixCSR,I,J,V,m,n)
julia> CSR= sparsecsr(I,J,V,m,n)

```

## Installation

**SparseMatricesCSR** itself is installed when you add and use it into another project.

To include into your project from Julia REPL, use the following commands:

```
pkg> add SparseMatricesCSR
julia> using SparseMatricesCSR
```

