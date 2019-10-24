
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
getindex(A::SymSparseMatrixCSR, x::Integer, y::Integer) = getindex(A.lowertrian,max(x,y),min(x,y))

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
    SymSparseMatrixCSR(sparse(J,I,V,m,n,combine))
end


"""
    function push_coo!(::Type{SymSparseMatrixCSR},I,J,V,ik,jk,vk) 

Inserts entries in COO vectors for further building a SymSparseMatrixCSR.
It stores only the upper triangle, ignoring entries with (ik>jk) coordinates.
"""
function push_coo!(::Type{SymSparseMatrixCSR},I::Vector,J::Vector,V::Vector,ik::Integer,jk::Integer,vk::Number) 
    (ik>jk) && return
    (push!(I, ik), push!(J, jk), push!(V, vk))
end

"""
    function finalize_coo!(::Type{SymSparseMatrixCSR},I,J,V,m,n) 

Check and insert diagonal entries in COO vectors if needed.
"""
function finalize_coo!(T::Type{SymSparseMatrixCSR},I::Vector,J::Vector,V::Vector, m::Integer, n::Integer)
    m == n || throw(DimensionMismatch("matrix is not square: dimensions are ($m, $n)"))
    touched = zeros(Bool,m)
    for k in 1:length(I)
        Ik = I[k]
        Jk = J[k]
        if Ik == Jk
            touched[Ik] = true
        end
    end
    for k in 1:m
        if ! touched[k]
            push_coo!(T,I,J,V,k,k,zero(eltype(V)))
        end
    end
end


