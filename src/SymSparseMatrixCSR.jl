
"""
    struct SymSparseMatrixCSR{T,Ti<:Integer} <: AbstractSparseMatrix{T,Ti}

Matrix type for storing symmetric sparse matrices in the
Compressed Sparse Row format. The standard way
of constructing SparseMatrixCSR is through the 
[`symsparsecsr`](@ref) function.
"""
struct SymSparseMatrixCSR{Bi,T,Ti<:Integer} <: AbstractSparseMatrix{T,Ti}
    uppertrian :: SparseMatrixCSR{Bi,T,Ti}
end


show(io::IO, A::SymSparseMatrixCSR) = show(io, A.uppertrian)
size(A::SymSparseMatrixCSR) = size(A.uppertrian)
getindex(A::SymSparseMatrixCSR, x::Integer, y::Integer) = getindex(A.uppertrian,min(x,y),max(x,y))

"""
    nnz(S::SymSparseMatrixCSR)

Returns the number of stored (filled) elements in a sparse array.
"""
nnz(S::SymSparseMatrixCSR) = nnz(S.uppertrian)

"""
    count(pred, S::SymSparseMatrixCSR) -> Integer

Count the number of elements in S for which predicate pred returns true. 
"""
count(pred, S::SymSparseMatrixCSR) = count(pred, S.uppertrian)

"""
    nonzeros(S::SymSparseMatrixCSR)

Return a vector of the structural nonzero values in sparse array S. 
This includes zeros that are explicitly stored in the sparse array. 
The returned vector points directly to the internal nonzero storage of A, 
and any modifications to the returned vector will mutate A as well.
"""
nonzeros(S::SymSparseMatrixCSR) = nonzeros(S.uppertrian)

"""
    nzrange(S::SymSparseMatrixCSR, row::Integer)

Return the range of indices to the structural nonzero values of a 
sparse matrix row. 
"""
nzrange(S::SymSparseMatrixCSR, row::Integer) = nzrange(S.uppertrian, row)

"""
    findnz(S::SymSparseMatrixCSR)

Return a tuple (I, J, V) where I and J are the row and column indices 
of the stored ("structurally non-zero") values in sparse matrix A, 
and V is a vector of the values.
"""
findnz(S::SymSparseMatrixCSR) = findnz(S.uppertrian)

"""
    rowvals(S::SymSparseMatrixCSR)

Return an error. 
CSR sparse matrices does not contain raw row values.
It contains row pointers instead that can be accessed
by using [`nzrange`](@ref).
"""

rowvals(S::SymSparseMatrixCSR) = rowvals(S.uppertrian)

"""
    colvals(S::SparseMatrixCSR)

Return a vector of the col indices of S. 
Any modifications to the returned vector will mutate S as well. 
Providing access to how the co,l indices are stored internally 
can be useful in conjunction with iterating over structural 
nonzero values. See also [`nonzeros`](@ref) and [`nzrange`](@ref).
"""

colvals(S::SymSparseMatrixCSR) = colvals(S.uppertrian)


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
function symsparsecsr(T::Type{SymSparseMatrixCSR{Bi}},I,J,V,args...)  where {Bi}
    m = length(args)>0 ? args[1] : isempty(I) ? 0 : Int(maximum(I))
    n = length(args)>1 ? args[2] : m
    c = length(args)>2 ? args[3] : +
    symsparsecsr(T,I,J,V,m,n,c)
end

function symsparsecsr(::Type{SymSparseMatrixCSR{Bi}},I,J,V,m,n,combine) where {Bi}
    m == n || throw(DimensionMismatch("matrix is not square: dimensions are ($m, $n)"))
    SymSparseMatrixCSR(sparsecsr(SparseMatrixCSR{Bi},I,J,V,m,n,combine))
end

symsparsecsr(I,J,V,args...) =
    symsparsecsr(SymSparseMatrixCSR{1},I,J,V,args...)


"""
    function push_coo!(::Type{SymSparseMatrixCSR},I,J,V,ik,jk,vk) 

Inserts entries in COO vectors for further building a SymSparseMatrixCSR.
It stores only the upper triangle, ignoring entries with (ik>jk) coordinates.
"""
function push_coo!(::Type{SymSparseMatrixCSR{Bi}},
        I::Vector,J::Vector,V::Vector,ik::Integer,jk::Integer,vk::Number) where {Bi}
    (ik>jk) && return
    (push!(I, ik), push!(J, jk), push!(V, vk))
end
push_coo!(::Type{SymSparseMatrixCSR}, args...) = push_coo!(SymSparseMatrixCSR{1}, args...) 

"""
    function finalize_coo!(::Type{SymSparseMatrixCSR},I,J,V,m,n) 

Check and insert diagonal entries in COO vectors if needed.
"""
function finalize_coo!(T::Type{SymSparseMatrixCSR{Bi}},
        I::Vector,J::Vector,V::Vector, m::Integer, n::Integer) where {Bi}
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

finalize_coo!(::Type{SymSparseMatrixCSR}, args...) = finalize_coo!(SymSparseMatrixCSR{1}, args...) 


"""
    function mul!(y::AbstractVector,A::SymSparseMatrixCSR,v::AbstractVector{T}) where {T}

Calculates the matrix-vector product ``Av`` and stores the result in `y`,
overwriting the existing value of `y`. 
"""
function mul!(y::AbstractVector,A::SymSparseMatrixCSR,v::AbstractVector{T}) where {T}
    A.uppertrian.n == size(v, 1) || throw(DimensionMismatch())
    A.uppertrian.m == size(y, 1) || throw(DimensionMismatch())

    y .= zero(T)
    for row = 1:size(y, 1)
        @inbounds for nz in nzrange(A,row)
            col = A.uppertrian.colval[nz]-A.uppertrian.offset
            y[row] += A.uppertrian.nzval[nz]*v[col]
            row != col && (y[col] += A.uppertrian.nzval[nz]*v[row])
        end
    end
    return y
end


