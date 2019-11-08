
"""
    struct SparseMatrixCSR{Idx,T,Ti<:Integer} <: AbstractSparseMatrix{T,Ti}

Matrix type for storing Idx-based sparse matrices 
in the Compressed Sparse Row format. The standard 
way of constructing SparseMatrixCSR is through the 
[`sparsecsr`](@ref) function.
"""
struct SparseMatrixCSR{Idx,T,Ti} <: AbstractSparseMatrix{T,Ti}
    transpose :: SparseMatrixCSC{T,Ti}
    function SparseMatrixCSR{Idx}(transpose::SparseMatrixCSC{T,Ti}) where {Idx,T,Ti}
        transpose.colptr .+= Idx-1
        transpose.rowval .+= Idx-1
        new{Idx,T,Ti}(transpose)
    end
end


SparseMatrixCSR(transpose::SparseMatrixCSC) = SparseMatrixCSR{1}(transpose)
show(io::IO, A::SparseMatrixCSR) = show(io, A.transpose)
size(A::SparseMatrixCSR) = (A.transpose.n, A.transpose.m)


function getindex(A::SparseMatrixCSR{Idx,T,Ti}, i0::Integer, i1::Integer) where {Idx,T,Ti}
    if !(Idx <= i0 <= A.transpose.n && Idx <= i1 <= A.transpose.m); throw(BoundsError()); end
    r1 = Int(A.transpose.colptr[i0]+1-Idx)
    r2 = Int(A.transpose.colptr[i0+1]-Idx)
    (r1 > r2) && return zero(T)
    i1_= i1-1+Idx
    r1 = searchsortedfirst(A.transpose.rowval, i1_, r1, r2, Base.Order.Forward)
    ((r1 > r2) || (A.transpose.rowval[r1] != i1_)) ? zero(T) : A.transpose.nzval[r1]
end

"""
    nnz(S::SparseMatrixCSR)

Returns the number of stored (filled) elements in a sparse array.
"""
nnz(S::SparseMatrixCSR{Idx}) where {Idx} = nnz(S.transpose)-Idx+1

"""
    count(pred, S::SparseMatrixCSR) -> Integer

Count the number of elements in S for which predicate pred returns true. 
"""
count(pred, S::SparseMatrixCSR) = count(pred, nonzeros(S))

"""
    nonzeros(S::SparseMatrixCSR)

Return a vector of the structural nonzero values in sparse array S. 
This includes zeros that are explicitly stored in the sparse array. 
The returned vector points directly to the internal nonzero storage of A, 
and any modifications to the returned vector will mutate A as well.
"""
nonzeros(S::SparseMatrixCSR) = S.transpose.nzval

"""
    nzrange(S::SparseMatrixCSR, row::Integer

Return the range of indices to the structural nonzero values of a 
sparse matrix row. 
"""
nzrange(S::SparseMatrixCSR{Idx}, row::Integer) where {Idx} = nzrange(S.transpose, row+1-Idx).+(1-Idx)

"""
    findnz(S::SparseMatrixCSR)

Return a tuple (I, J, V) where I and J are the row and column indices 
of the stored ("structurally non-zero") values in sparse matrix A, 
and V is a vector of the values.
"""
function findnz(S::SparseMatrixCSR{Idx,Tv,Ti}) where {Idx,Tv,Ti}
    numnz = nnz(S)
    I = Vector{Ti}(undef, numnz)
    J = Vector{Ti}(undef, numnz)
    V = Vector{Tv}(undef, numnz)

    count = 1
    offset = 1-Idx
    @inbounds for col = Idx : S.transpose.n-offset, k = nzrange(S,col)
        I[count] = S.transpose.rowval[k]+offset
        J[count] = col+offset
        V[count] = S.transpose.nzval[k]
        count += 1
    end

    return (I, J, V)
end

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
sparsecsr(::Type{SparseMatrixCSR}, I,J,args...) = SparseMatrixCSR(sparse(J,I,args...))
sparsecsr(::Type{SparseMatrixCSR{idx}}, I,J,args...) where {idx} = SparseMatrixCSR{idx}(sparse(J,I,args...))


"""
    function push_coo!(::Type{SparseMatrixCSR},I,J,V,ik,jk,vk) 

Inserts entries in COO vectors for further building a SparseMatrixCSR.
"""
function push_coo!(::Type{SparseMatrixCSR{Idx}},
    I::Vector,J::Vector,V::Vector,ik::Integer,jk::Integer,vk::Number) where {Idx}
    (push!(I, ik), push!(J, jk), push!(V, vk))
end
push_coo!(SparseMatrixCSR, args...) = push_coo!(SparseMatrixCSR{1}, args...) 

"""
    function finalize_coo!(::Type{SparseMatrixCSR},I,J,V,m,n) 

Check and insert diagonal entries in COO vectors if needed.
"""
function finalize_coo!(::Type{SparseMatrixCSR{Idx}},
    I::Vector,J::Vector,V::Vector,m::Integer,n::Integer) where {Idx}
end

finalize_coo!(SparseMatrixCSR, args...) = finalize_coo!(SparseMatrixCSR{1}, args...) 




