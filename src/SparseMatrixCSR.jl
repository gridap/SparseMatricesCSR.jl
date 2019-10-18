
"""
    struct SparseMatrixCSR{T,Ti<:Integer} <: AbstractSparseMatrix{T,Ti}

Matrix type for storing sparse matrices in the
Compressed Sparse Row format. The standard way
of constructing SparseMatrixCSR is through the 
[`sparsecsr`](@ref) function.
"""
struct SparseMatrixCSR{T,Ti<:Integer} <: AbstractSparseMatrix{T,Ti}
    transpose :: SparseMatrixCSC{T,Ti}
end


show(io::IO, A::SparseMatrixCSR) = show(io, A.transpose)
size(A::SparseMatrixCSR) = (A.transpose.n, A.transpose.m)
getindex(A::SparseMatrixCSR, x::Integer, y::Integer) = getindex(A.transpose,y,x)

"""
    nnz(S::SparseMatrixCSR)

Returns the number of stored (filled) elements in a sparse array.
"""
nnz(S::SparseMatrixCSR) = nnz(S.transpose)

"""
    count(pred, S::SparseMatrixCSR) -> Integer

Count the number of elements in S for which predicate pred returns true. 
"""
count(pred, S::SparseMatrixCSR) = count(pred, S.transpose)

"""
    nonzeros(S::SparseMatrixCSR)

Return a vector of the structural nonzero values in sparse array S. 
This includes zeros that are explicitly stored in the sparse array. 
The returned vector points directly to the internal nonzero storage of A, 
and any modifications to the returned vector will mutate A as well.
"""
nonzeros(S::SparseMatrixCSR) = S.transpose.nzval

"""
    nzrange(S::SparseMatrixCSR, row::Integer)

Return the range of indices to the structural nonzero values of a 
sparse matrix row. 
"""
nzrange(S::SparseMatrixCSR, row::Integer) = nzrange(S.transpose, row)

"""
    findnz(S::SparseMatrixCSR)

Return a tuple (I, J, V) where I and J are the row and column indices 
of the stored ("structurally non-zero") values in sparse matrix A, 
and V is a vector of the values.
"""
findnz(S::SparseMatrixCSR) = findnz(S.transpose)

"""
    sparsecsr(I, J, V, [m, n, combine])

Create a sparse matrix S of dimensions m x n such that S[I[k], J[k]] = V[k]. 
The combine function is used to combine duplicates. 
If m and n are not specified, they are set to 
maximum(I) and maximum(J) respectively. 
If the combine function is not supplied, combine defaults to + 
unless the elements of V are Booleans in which case combine defaults to |. 
All elements of I must satisfy 1 <= I[k] <= m, 
and all elements of J must satisfy 1 <= J[k] <= n. 
Numerical zeros in (I, J, V) are retained as structural nonzeros; 
to drop numerical zeros, use dropzeros!.
"""
sparsecsr(I,J,args...) = SparseMatrixCSR(sparse(J,I,args...))

"""
    function push_coo!(::Type{SparseMatrixCSR},I,J,V,ik,jk,vk) 

Inserts entries in COO vectors for further building a SparseMatrixCSR.
"""
function push_coo!(::Type{SparseMatrixCSR},I,J,V,ik,jk,vk) 
    (push!(I, jk), push!(J, ik), push!(V, vk))
end

