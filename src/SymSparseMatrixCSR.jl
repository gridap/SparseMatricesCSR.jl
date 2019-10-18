
"""
    struct SymSparseMatrixCSR{T,Ti<:Integer} <: AbstractSparseMatrix{T,Ti}

Matrix type for storing symmetric sparse matrices in the
Compressed Sparse Row format. The standard way
of constructing SparseMatrixCSR is through the 
[`symsparsecsr`](@ref) function.
"""
struct SymSparseMatrixCSR{T,Ti<:Integer} <: AbstractSparseMatrix{T,Ti}
    lowertrian :: SparseMatrixCSC{T,Ti}
end


show(io::IO, A::SymSparseMatrixCSR) = show(io, A.lowertrian)
size(A::SymSparseMatrixCSR) = size(A.lowertrian)
getindex(A::SymSparseMatrixCSR, x::Integer, y::Integer) = getindex(A.lowertrian,min(x,y),max(x,y))

"""
    nnz(S::SymSparseMatrixCSR)

Returns the number of stored (filled) elements in a sparse array.
"""
nnz(S::SymSparseMatrixCSR) = nnz(S.lowertrian)

"""
    count(pred, S::SymSparseMatrixCSR) -> Integer

Count the number of elements in S for which predicate pred returns true. 
"""
count(pred, S::SymSparseMatrixCSR) = count(pred, S.lowertrian)

"""
    nonzeros(S::SymSparseMatrixCSR)

Return a vector of the structural nonzero values in sparse array S. 
This includes zeros that are explicitly stored in the sparse array. 
The returned vector points directly to the internal nonzero storage of A, 
and any modifications to the returned vector will mutate A as well.
"""
nonzeros(S::SymSparseMatrixCSR) = S.lowertrian.nzval

"""
    nzrange(S::SymSparseMatrixCSR, row::Integer)

Return the range of indices to the structural nonzero values of a 
sparse matrix row. 
"""
nzrange(S::SymSparseMatrixCSR, row::Integer) = nzrange(S.lowertrian, row)

"""
    findnz(S::SymSparseMatrixCSR)

Return a tuple (I, J, V) where I and J are the row and column indices 
of the stored ("structurally non-zero") values in sparse matrix A, 
and V is a vector of the values.
"""
findnz(S::SymSparseMatrixCSR) = findnz(S.lowertrian)

"""
    symsparsecsr(I, J, V, [m, n, combine])

Create a Symmetric sparse matrix S of dimensions m x n
 such that S[I[k], J[k]] = V[k], and m=n. 
The combine function is used to combine duplicates. 
If m and n are not specified, they are set to 
maximum(I) and maximum(J) respectively. 
If the combine function is not supplied, combine defaults to +.
All elements of I must satisfy 1 <= I[k] <= m, 
and all elements of J must satisfy 1 <= J[k] <= n. 
Numerical zeros in (I, J, V) are retained as structural nonzeros; 
to drop numerical zeros, use dropzeros!.
"""
function symsparsecsr(I,J,V,args...) 
    m = length(args)>0 ? args[1] : isempty(I) ? 0 : Int(maximum(I))
    n = length(args)>1 ? args[2] : m
    c = length(args)>2 ? args[3] : +
    symsparsecsr(I,J,V,m,n,c)
end

function symsparsecsr(I,J,V,m,n,combine) 
    m == n || throw(DimensionMismatch("matrix is not square: dimensions are ($m, $n)"))
    indices = [ i for i in 1:length(I) if !(I[i]<J[i])]
    # Explicitly store diagonal zeros if needed
    SymSparseMatrixCSR(sparse(vcat(J[indices],1:m),vcat(I[indices],1:m),vcat(V[indices],zeros(eltype(V),m)),m,m,combine))
end

