
"""
    struct SymSparseMatrixCSR{Bi,T,Ti<:Integer} <: AbstractSparseMatrix{T,Ti}
      uppertrian :: SparseMatrixCSR{Bi,T,Ti}
    end

Matrix type for storing symmetric sparse matrices in the
Compressed Sparse Row format with `Bi`-based indexing (typically 0 or 1).
Only the upper triangle is stored (including the non zero diagonal entries),
which is represented by a `SparseMatrixCSR`.
The standard way of constructing a `SymSparseMatrixCSR` is through the
[`symsparsecsr`](@ref) function.
"""
struct SymSparseMatrixCSR{Bi,T,Ti<:Integer} <: AbstractSparseMatrix{T,Ti}
    uppertrian :: SparseMatrixCSR{Bi,T,Ti}
end

"""
    symsparsecsr(args...;symmetrize::Bool=false)
    symsparsecsr(::Val{Bi},args...;symmetrize::Bool=false) where Bi

Create  a `SymSparseMatrixCSR` with `Bi`-based indexing (1 by default)
from the same `args...` as one constructs a `SparseMatrixCSC`
with the [`sparse`](@extref SparseArrays.sparse) function. If `symmetrize == false` (the default)
the given arguments should only describe the upper triangle
of the matrix (including non zero diagonal values). If `symmetrize == true`
a non symmetric input is accepted and it will be symmetrized in-place
(i.e., changing the input arguments).
"""
function symsparsecsr(::Val{Bi},I,J,V,args...;symmetrize::Bool=false) where Bi
  if symmetrize
    Tv = eltype(V)
    α = Tv(0.5)
    for k in 1:length(I)
      r = I[k]
      c = J[k]
      if r > c
        I[k] = c
        J[k] = r
      end
      if r != c
        V[k] = α*V[k]
      end
    end
  end
  SymSparseMatrixCSR(sparsecsr(Val(Bi),I,J,V,args...))
end
symsparsecsr(args...;kwargs...) = symsparsecsr(Val(1),args...;kwargs...)

size(A::SymSparseMatrixCSR) = size(A.uppertrian)
IndexStyle(::Type{<:SymSparseMatrixCSR}) = IndexCartesian()
function getindex(A::SymSparseMatrixCSR, x::Integer, y::Integer)
  getindex(A.uppertrian,min(x,y),max(x,y))
end
function setindex!(A::SymSparseMatrixCSR, v, x::Integer, y::Integer)
  setindex!(A.uppertrian,v,min(x,y),max(x,y))
end

getrowptr(S::SymSparseMatrixCSR) = getrowptr(S.uppertrian)
getnzval(S::SymSparseMatrixCSR) = getnzval(S.uppertrian)
getcolval(S::SymSparseMatrixCSR) = getcolval(S.uppertrian)

"""
    getBi(S::SymSparseMatrixCSR{Bi}) where {Bi}

Return `Bi`.
"""
getBi(S::SymSparseMatrixCSR{Bi}) where {Bi} = Bi

"""
    getoffset(S::SymSparseMatrixCSR{Bi}) where {Bi}

Return `1-Bi`. Useful to convert from 1-based to `Bi`-based indexing
(by subtracting the offset).
"""
getoffset(S::SymSparseMatrixCSR{Bi}) where Bi = getoffset(Bi)

"""
    issparse(S::SymSparseMatrixCSR)

Returns `true`.
"""
issparse(S::SymSparseMatrixCSR) = true

"""
    nnz(S::SymSparseMatrixCSR)

Returns the number of stored elements in a sparse array,
which correspond to the nonzero entries in the upper triangle and diagonal.
"""
nnz(S::SymSparseMatrixCSR) = nnz(S.uppertrian)

"""
    nonzeros(S::SymSparseMatrixCSR)

Return a vector (1-based) of the structural nonzero values in sparse array S.
This includes zeros that are explicitly stored in the sparse array,
which correspond to the nonzero entries in the upper triangle and diagonal.
The returned vector points directly to the internal nonzero storage of S,
and any modifications to the returned vector will mutate S as well.
"""
nonzeros(S::SymSparseMatrixCSR) = nonzeros(S.uppertrian)

"""
    colvals(S::SparseMatrixCSR)

Return a vector of the col indices of `S`. The stored values are indexes to arrays
with `Bi`-based indexing, but the `colvals(S)` array itself is a standard 1-based
Julia `Vector`.
Any modifications to the returned vector will mutate S as well.
Providing access to how the col indices are stored internally
can be useful in conjunction with iterating over structural
nonzero values. See also [`nonzeros`](@ref) and [`nzrange`](@ref).
"""
colvals(S::SymSparseMatrixCSR) = colvals(S.uppertrian)

"""
    nzrange(S::SymSparseMatrixCSR, row::Integer)

Return the range of indices to the structural nonzero values of a
sparse matrix row section being in the diagonal or upper triangle.
The returned range of indices is always 1-based even for `Bi != 1`.
"""
nzrange(S::SymSparseMatrixCSR, row::Integer) = nzrange(S.uppertrian, row)

"""
    findnz(S::SymSparseMatrixCSR)

Return a tuple `(I, J, V)` where `I` and `J` are the row and column 1-based indices
of the stored ("structurally non-zero in diagonal + upper trianle") values in sparse matrix A,
and V is a vector of the values. The returned vectors are newly allocated
and are unrelated to the internal storage of matrix `S`.
"""
findnz(S::SymSparseMatrixCSR) = findnz(S.uppertrian)

"""
    count(pred, S::SymSparseMatrixCSR)
    count(S::SymSparseMatrixCSR)

Count the number of elements in `nonzeros(S)` for which predicate `pred` returns `true`.
If  `pred` not given, it counts the number of `true` values.
"""
count(pred, S::SymSparseMatrixCSR) = count(pred, S.uppertrian)
count(S::SymSparseMatrixCSR) = count(i->true, S)

function LinearAlgebra.fillstored!(a::SymSparseMatrixCSR,v)
  LinearAlgebra.fillstored!(a.uppertrian,v)
  a
end

function LinearAlgebra.rmul!(a::SymSparseMatrixCSR,v::Number)
  LinearAlgebra.rmul!(a.uppertrian,v)
  a
end

function mul!(y::AbstractVector,A::SymSparseMatrixCSR,v::AbstractVector, α::Number, β::Number)
  A.uppertrian.n == size(v, 1) || throw(DimensionMismatch())
  A.uppertrian.m == size(y, 1) || throw(DimensionMismatch())
  if β != 1
    β != 0 ? rmul!(y, β) : fill!(y, zero(eltype(y)))
  end
  o = getoffset(A)
  for row = 1:size(y, 1)
    @inbounds for nz in nzrange(A,row)
      col = A.uppertrian.colval[nz]+o
      y[row] += A.uppertrian.nzval[nz]*v[col]*α
      row != col && (y[col] += A.uppertrian.nzval[nz]*v[row]*α)
    end
  end
  return y
end

function mul!(y::AbstractVector,A::SymSparseMatrixCSR,v::AbstractVector)
  A.uppertrian.n == size(v, 1) || throw(DimensionMismatch())
  A.uppertrian.m == size(y, 1) || throw(DimensionMismatch())

  fill!(y,zero(eltype(y)))
  o = getoffset(A)
  for row = 1:size(y, 1)
    @inbounds for nz in nzrange(A,row)
      col = A.uppertrian.colval[nz]+o
      y[row] += A.uppertrian.nzval[nz]*v[col]
      row != col && (y[col] += A.uppertrian.nzval[nz]*v[row])
    end
  end
  return y
end

*(A::SymSparseMatrixCSR, v::Vector) = (y = similar(v,size(A,1));mul!(y,A,v))

function show(io::IO, ::MIME"text/plain", S::SymSparseMatrixCSR)
    xnnz = nnz(S)
    print(io, S.uppertrian.m, "×", S.uppertrian.n, " ",
              typeof(S), " with ", xnnz, " stored ",
              xnnz == 1 ? "entry" : "entries")
    if xnnz != 0
        print(io, ":")
        show(IOContext(io, :typeinfo => eltype(S)), S)
    end
end
show(io::IO, S::SymSparseMatrixCSR) = show(io, S.uppertrian)

Base.convert(::Type{T},a::T) where T<:SymSparseMatrixCSR = a
function Base.convert(
  ::Type{SymSparseMatrixCSR{Bi,Tv,Ti}},a::SymSparseMatrixCSR) where {Bi,Tv,Ti}
  utrian = convert(SparseMatrixCSR{Bi,Tv,Ti},a.uppertrian)
  SymSparseMatrixCSR(utrian)
end

function Base.copy(a::SymSparseMatrixCSR{Bi,T,Ti}) where {Bi,T,Ti}
   SymSparseMatrixCSR{Bi,T,Ti}(copy(a.uppertrian))
end
