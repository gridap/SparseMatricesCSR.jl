
"""
    struct SparseMatrixCSR{Bi,Tv,Ti<:Integer} <: AbstractSparseMatrix{Tv,Ti}

Matrix type for storing Bi-based sparse matrices 
in the Compressed Sparse Row format. The standard 
way of constructing SparseMatrixCSR is through the 
[`sparsecsr`](@ref) function.
"""
struct SparseMatrixCSR{Bi,Tv,Ti} <: AbstractSparseMatrix{Tv,Ti}
    offset :: Int             # Bi-index offset
    m      :: Int             # Number of columns
    n      :: Int             # Number of rows
    rowptr :: Vector{Ti}      # Row i is in rowptr[i]:(rowptr[i+1]-1)
    colval :: Vector{Ti}      # Col indices of stored values
    nzval  :: Vector{Tv}      # Stored values, typically nonzeros

    function SparseMatrixCSR{Bi,Tv,Ti}(m::Integer, n::Integer, rowptr::Vector{Ti}, colval::Vector{Ti},
                                    nzval::Vector{Tv}) where {Bi,Tv,Ti<:Integer}
        @noinline throwsz(str, lbl, k) =
            throw(ArgumentError("number of $str ($lbl) must be ≥ 0, got $k"))
        m < 0 && throwsz("rows", 'm', m)
        n < 0 && throwsz("columns", 'n', n)
        offset = Bi-1
        rowptr .+= offset
        colval .+= offset
        new{Bi,Tv,Ti}(Int(offset), Int(m), Int(n), rowptr, colval, nzval)
    end
end


SparseMatrixCSR(transpose::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti} = 
        SparseMatrixCSR{1,Tv,Ti}(transpose.n, transpose.m, transpose.colptr, transpose.rowval, transpose.nzval)
SparseMatrixCSR{Bi}(transpose::SparseMatrixCSC{Tv,Ti}) where {Bi,Tv,Ti} = 
        SparseMatrixCSR{Bi,Tv,Ti}(transpose.n, transpose.m, transpose.colptr, transpose.rowval, transpose.nzval)
show(io::IO, S::SparseMatrixCSR) = show(io, S)
size(S::SparseMatrixCSR) = (S.m, S.n)


function getindex(S::SparseMatrixCSR{Bi,Tv,Ti}, i0::Integer, i1::Integer) where {Bi,Tv,Ti}
    if !(Bi <= i0 <= S.m && Bi <= i1 <= S.n); throw(BoundsError()); end
    r1 = Int(S.rowptr[i0]-S.offset)
    r2 = Int(S.rowptr[i0+1]-Bi)
    (r1 > r2) && return zero(Tv)
    i1_= i1+S.offset
    r1 = searchsortedfirst(S.colval, i1_, r1, r2, Base.Order.Forward)
    ((r1 > r2) || (S.colval[r1] != i1_)) ? zero(Tv) : S.nzval[r1]
end

"""
    nnz(S::SparseMatrixCSR{Bi}) where {Bi}

Returns the number of stored (filled) elements in a sparse array.
"""
nnz(S::SparseMatrixCSR{Bi}) where {Bi} = S.rowptr[S.m+1]-Bi

"""
    count(pred, S::SparseMatrixCSR) -> Integer

Count the number of elements in S for which predicate pred returns true. 
"""
count(pred, S::SparseMatrixCSR) = count(pred, nonzeros(S))

"""
    nonzeros(S::SparseMatrixCSR)

Return a vector of the structural nonzero values in sparse array S. 
This includes zeros that are explicitly stored in the sparse array. 
The returned vector points directly to the internal nonzero storage of S, 
and any modifications to the returned vector will mutate S as well.
"""
nonzeros(S::SparseMatrixCSR) = S.nzval

"""
    nzrange(S::SparseMatrixCSR{Bi}, row::Integer) where {Bi}

Return the range of indices to the structural nonzero values of a 
sparse matrix row. 
"""
nzrange(S::SparseMatrixCSR{Bi}, row::Integer) where {Bi} = S.rowptr[row]-S.offset:S.rowptr[row+1]-Bi

"""
    function findnz(S::SparseMatrixCSR{Bi,Tv,Ti})

Return a tuple (I, J, V) where I and J are the row and column indices 
of the stored ("structurally non-zero") values in sparse matrix A, 
and V is a vector of the values.
"""
function findnz(S::SparseMatrixCSR{Bi,Tv,Ti}) where {Bi,Tv,Ti}
    numnz = nnz(S)
    I = Vector{Ti}(undef, numnz)
    J = Vector{Ti}(undef, numnz)
    V = Vector{Tv}(undef, numnz)

    count = 1
    @inbounds for row = 1 : S.m, k = nzrange(S,row)
        I[count] = S.colval[k]-S.offset
        J[count] = row
        V[count] = S.nzval[k]
        count += 1
    end

    return (I, J, V)
end

"""
    rowvals(S::SparseMatrixCSR)

Return an error. 
CSR sparse matrices does not contain raw row values.
It contains row pointers instead that can be accessed
by using [`nzrange`](@ref).
"""
rowvals(S::SparseMatrixCSR) = error("CSR sparse matrix does not contain raw row values")

"""
    colvals(S::SparseMatrixCSR)

Return a vector of the col indices of S. 
Any modifications to the returned vector will mutate S as well. 
Providing access to how the co,l indices are stored internally 
can be useful in conjunction with iterating over structural 
nonzero values. See also [`nonzeros`](@ref) and [`nzrange`](@ref).
"""

colvals(S::SparseMatrixCSR) = S.colval

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
sparsecsr(I,J,args...) = 
        SparseMatrixCSR(sparse(J,I,args...))
sparsecsr(::Type{SparseMatrixCSR},I,J,args...) = 
        sparsecsr(I,J,args...)
sparsecsr(::Type{SparseMatrixCSR{Bi}},I,J,args...) where {Bi} = 
        SparseMatrixCSR{Bi}(sparse(J,I,args...))
sparsecsr(::Type{SparseMatrixCSR{Bi}},I,J,V,m,n,args...) where {Bi} = 
        SparseMatrixCSR{Bi}(sparse(J,I,V,n,m,args...))


"""
    function push_coo!(::Type{SparseMatrixCSR},I,J,V,ik,jk,vk) 

Inserts entries in COO vectors for further building a SparseMatrixCSR.
"""
function push_coo!(::Type{SparseMatrixCSR{Bi}},
    I::Vector,J::Vector,V::Vector,ik::Integer,jk::Integer,vk::Number) where {Bi}
    (push!(I, ik), push!(J, jk), push!(V, vk))
end
push_coo!(::Type{SparseMatrixCSR}, args...) = push_coo!(SparseMatrixCSR{1}, args...) 

"""
    function finalize_coo!(::Type{SparseMatrixCSR},I,J,V,m,n) 

Check and insert diagonal entries in COO vectors if needed.
"""
function finalize_coo!(::Type{SparseMatrixCSR{Bi}},
    I::Vector,J::Vector,V::Vector,m::Integer,n::Integer) where {Bi}
end
finalize_coo!(::Type{SparseMatrixCSR}, args...) = finalize_coo!(SparseMatrixCSR{1}, args...) 



"""
    function mul!(y::AbstractVector,A::SparseMatrixCSR,v::AbstractVector{T}) where {T}

Calculates the matrix-vector product ``Av`` and stores the result in `y`,
overwriting the existing value of `y`. 
"""
function mul!(y::AbstractVector,A::SparseMatrixCSR,v::AbstractVector{T}) where {T}
    A.n == size(v, 1) || throw(DimensionMismatch())
    A.m == size(y, 1) || throw(DimensionMismatch())

    y .= zero(T)
    for row = 1:size(y, 1)
        @inbounds for nz in nzrange(A,row)
            col = A.colval[nz]-A.offset
            y[row] += A.nzval[nz]*v[col]
        end
    end
    return y
end

