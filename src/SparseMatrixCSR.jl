
"""
    SparseMatrixCSR{Bi,Tv,Ti<:Integer} <: AbstractSparseMatrix{Tv,Ti}

Matrix type for storing sparse matrices
in the Compressed Sparse Row format with `Bi`-based indexing (typically 0 or 1).
The standard way of constructing `SparseMatrixCSR` is through the
[`sparsecsr`](@ref) function.

# Properties

- `m::Int` Number of columns
- `n::Int` Number of rows
- `rowptr::Vector{Ti}` Row `i` (1-based) in `nzval` vector at indices (1-based) `(rowptr[i]+(1-Bi)):(rowptr[i+1]+(1-Bi)-1)`
- `colval::Vector{Ti}` Col indices (`Bi`-based) of stored values
- `nzval::Vector{Tv}`  Stored values, typically non zeros

# Inner constructor

    SparseMatrixCSR{Bi}(
      m::Integer,
      n::Integer,
      rowptr::Vector{Ti},
      colval::Vector{Ti},
      nzval::Vector{Tv}) where {Bi,Tv,Ti<:Integer}

"""
struct SparseMatrixCSR{Bi,Tv,Ti} <: AbstractSparseMatrix{Tv,Ti}
  m::Int
  n::Int
  rowptr::Vector{Ti}
  colval::Vector{Ti}
  nzval::Vector{Tv}
  function SparseMatrixCSR{Bi}(m::Integer, n::Integer, rowptr::Vector{Ti}, colval::Vector{Ti},
    nzval::Vector{Tv}) where {Bi,Tv,Ti<:Integer}
    @noinline throwsz(str, k, r) =
      throw(ArgumentError("$str must be $r, got $k"))
    m < 0 && throwsz("m",">=0",m)
    n < 0 && throwsz("n",">=0", n)
    length(rowptr) != m+1 && throwsz("lengh(rowptr)",m+1,length(rowptr))
    rowptr[end]-Bi != length(nzval) && throwsz("rowptr[end]-Bi",rowptr[m+1]-Bi, length(nzval))
    new{Bi,Tv,Ti}(Int(m), Int(n), rowptr, colval, nzval)
  end
end

"""
    SparseMatrixCSR(a::Transpose{Tv,<:SparseMatrixCSC} where Tv)

Build a 1-based `SparseMatrixCSR` from the lazy transpose of a `SparseMatrixCSC`.
The resulting matrix takes ownership of the internal storage of input matrix.
Modifying the values of one, will mutate also the other.
"""
function SparseMatrixCSR(a::Transpose{Tv,<:SparseMatrixCSC} where Tv)
  at = a.parent
  SparseMatrixCSR{1}(size(a,1),size(a,2),at.colptr,rowvals(at),nonzeros(at))
end

"""
    SparseMatrixCSR(a::SparseMatrixCSC}

Build a 1-based `SparseMatrixCSR` from a `SparseMatrixCSC`. 
"""
SparseMatrixCSR(a::SparseMatrixCSC) = SparseMatrixCSR(transpose(sparse(transpose(a))))

"""
    SparseMatrixCSR(a::AbstractMatrix}

Build a 1-based `SparseMatrixCSR` from an `AbstractMatrix`. 
"""
SparseMatrixCSR(a::AbstractMatrix) = SparseMatrixCSR(sparse(a))

"""
    SparseMatrixCSR{Bi}(a::Transpose{Tv,<:SparseMatrixCSC} where Tv) where Bi

Build a `Bi`-based `SparseMatrixCSR` from the lazy transpose of a `SparseMatrixCSC`.
The resulting matrix takes ownership of the internal storage of input matrix
for any `Bi` and modifies the internal storage when `Bi != 1`.
The input matrix will become unusable in the latter case.
"""
function SparseMatrixCSR{Bi}(a::Transpose{Tv,<:SparseMatrixCSC} where Tv) where Bi
  at = a.parent
  at.colptr .-= getoffset(Bi)
  rowvals(at) .-= getoffset(Bi)
  SparseMatrixCSR{Bi}(size(a,1),size(a,2),at.colptr,rowvals(at),nonzeros(at))
end
SparseMatrixCSR{1}(a::Transpose{Tv,<:SparseMatrixCSC} where Tv) = SparseMatrixCSR(a)

"""
    sparsecsr(args...)
    sparsecsr(::Val{Bi},args...) where Bi

Create  a `SparseMatrixCSR` with `Bi`-based indexing (1 by default)
from the same `args...` as one constructs a `SparseMatrixCSC`
with the [`SparseArrays.sparse`](@extref) function.
"""
sparsecsr(A::AbstractMatrix{T}) where T = convert(SparseMatrixCSR{1,T,Int64}, A)
sparsecsr(I,J,V) = SparseMatrixCSR(transpose(sparse(J,I,V,dimlub(J),dimlub(I))))
sparsecsr(I,J,V,m,n) = SparseMatrixCSR(transpose(sparse(J,I,V,n,m)))
sparsecsr(I,J,V,m,n,combine) = SparseMatrixCSR(transpose(sparse(J,I,V,n,m,combine)))
sparsecsr(::Val{Bi},A::AbstractMatrix{T}) where {Bi,T} = convert(SparseMatrixCSR{Bi,T,Int64}, A)
sparsecsr(::Val{Bi},I,J,V) where Bi = SparseMatrixCSR{Bi}(transpose(sparse(J,I,V,dimlub(J),dimlub(I))))
sparsecsr(::Val{Bi},I,J,V,m,n) where Bi = SparseMatrixCSR{Bi}(transpose(sparse(J,I,V,n,m)))
sparsecsr(::Val{Bi},I,J,V,m,n,combine) where Bi = SparseMatrixCSR{Bi}(transpose(sparse(J,I,V,n,m,combine)))
dimlub(I) = isempty(I) ? 0 : Int(maximum(I))

Base.convert(::Type{T},a::T) where T<:SparseMatrixCSR = a
function Base.convert(
  ::Type{SparseMatrixCSR{Bi,Tv,Ti}},a::SparseMatrixCSR{Bi}) where {Bi,Tv,Ti}
  rowptr = convert(Vector{Ti},a.rowptr)
  colval = convert(Vector{Ti},a.colval)
  nzval = convert(Vector{Tv},a.nzval)
  SparseMatrixCSR{Bi}(a.m,a.n,rowptr,colval,nzval)
end
function Base.convert(
  ::Type{SparseMatrixCSR{Bi,Tv,Ti}},a::SparseMatrixCSR{Bj}) where {Bi,Tv,Ti,Bj}
  rowptr = Vector{Ti}(undef,length(a.rowptr))
  colval = Vector{Ti}(undef,length(a.colval))
  rowptr .= a.rowptr .+ (Bi-Bj)
  colval .= a.colval .+ (Bi-Bj)
  nzval = convert(Vector{Tv},a.nzval)
  SparseMatrixCSR{Bi}(a.m,a.n,rowptr,colval,nzval)
end
function Base.convert(
  ::Type{SparseMatrixCSR{Bi,Tv,Ti}},a::Transpose{Tu,<:SparseMatrixCSC} where Tu) where {Bi,Tv,Ti}
  convert(SparseMatrixCSR{Bi,Tv,Ti},SparseMatrixCSR{Bi}(a))
end
function Base.convert(
  ::Type{SparseMatrixCSR{Bi,Tv,Ti}},a::AbstractMatrix) where {Bi,Tv,Ti}
  at = convert(SparseMatrixCSC{Tv,Ti},transpose(a))
  convert(SparseMatrixCSR{Bi,Tv,Ti},transpose(at))
end

function Base.copy(a::SparseMatrixCSR{Bi}) where Bi
  SparseMatrixCSR{Bi}(a.m,a.n,copy(a.rowptr),copy(a.colval),copy(a.nzval))
end

_copy_and_increment(x) = copy(x) .+ 1

function LinearAlgebra.fillstored!(a::SparseMatrixCSR,v)
  fill!(a.nzval,v)
  a
end

function LinearAlgebra.rmul!(a::SparseMatrixCSR,v::Number)
  rmul!(a.nzval,v)
  a
end

function LinearAlgebra.lu(a::SparseMatrixCSR{0})
  rowptr = _copy_and_increment(a.rowptr)
  colval = _copy_and_increment(a.colval)
  transpose(lu(SparseMatrixCSC(a.m,a.n,rowptr,colval,a.nzval)))
end

function LinearAlgebra.lu(a::SparseMatrixCSR{1})
  transpose(lu(SparseMatrixCSC(a.m,a.n,a.rowptr,a.colval,a.nzval)))
end

if Base.USE_GPL_LIBS

const TransposeFact = isdefined(LinearAlgebra, :TransposeFactorization) ?
  LinearAlgebra.TransposeFactorization :
  Transpose

function LinearAlgebra.lu!(
    translu::TransposeFact{T,<:SuiteSparse.UMFPACK.UmfpackLU{T}},
    a::SparseMatrixCSR{1}) where {T}
  transpose(lu!(translu.parent,SparseMatrixCSC(a.m,a.n,a.rowptr,a.colval,a.nzval)))
end

function LinearAlgebra.lu!(
  translu::TransposeFact{T,<:SuiteSparse.UMFPACK.UmfpackLU{T}},
  a::SparseMatrixCSR{0}) where {T}
  rowptr = _copy_and_increment(a.rowptr)
  colval = _copy_and_increment(a.colval)
  transpose(lu!(translu.parent,SparseMatrixCSC(a.m,a.n,rowptr,colval,a.nzval)))
end

end # Base.USE_GPL_LIBS

size(S::SparseMatrixCSR) = (S.m, S.n)
IndexStyle(::Type{<:SparseMatrixCSR}) = IndexCartesian()
function getindex(A::SparseMatrixCSR{Bi,T}, i0::Integer, i1::Integer) where {Bi,T}
  if !(1 <= i0 <= size(A, 1) && 1 <= i1 <= size(A, 2)); throw(BoundsError()); end
  o = getoffset(A)
  r1 = Int(getrowptr(A)[i0]+o)
  r2 = Int(getrowptr(A)[i0+1]-Bi)
  (r1 > r2) && return zero(T)
  i1o = i1-o
  k = searchsortedfirst(colvals(A), i1o, r1, r2, Base.Order.Forward)
  ((k > r2) || (colvals(A)[k] != i1o)) ? zero(T) : nonzeros(A)[k]
end
function setindex!(A::SparseMatrixCSR{Bi,Tv,Ti}, _v, _i0::Integer, _i1::Integer) where {Bi,Tv,Ti}
  errmsg = "Trying to set an entry outside sparsity pattern"
  v  = convert(Tv,_v)
  i0 = convert(Ti,_i0)
  i1 = convert(Ti,_i1)
  if !(1 <= i0 <= size(A, 1) && 1 <= i1 <= size(A, 2)); throw(BoundsError()); end
  o = getoffset(A)
  r1 = Int(getrowptr(A)[i0]+o)
  r2 = Int(getrowptr(A)[i0+1]-Bi)
  (r1 > r2) && throw(ArgumentError(errmsg))
  i1o = i1-o
  k = searchsortedfirst(colvals(A), i1o, r1, r2, Base.Order.Forward)
  if ((k > r2) || (colvals(A)[k] != i1o))
    throw(ArgumentError(errmsg))
  end
  A.nzval[k]=v
end




getrowptr(S::SparseMatrixCSR) = S.rowptr
getnzval(S::SparseMatrixCSR) = S.nzval
getcolval(S::SparseMatrixCSR) = S.colval

"""
    getBi(S::SparseMatrixCSR{Bi}) where {Bi}

Return `Bi`.
"""
getBi(S::SparseMatrixCSR{Bi}) where {Bi} = Bi

"""
    getoffset(S::SparseMatrixCSR{Bi}) where {Bi}

Return `1-Bi`. Useful to convert from 1-based to `Bi`-based indexing
(by subtracting the offset).
"""
getoffset(S::SparseMatrixCSR{Bi}) where Bi = getoffset(Bi)
@inline getoffset(Bi::Integer) = 1-Bi

"""
    issparse(S::SparseMatrixCSR)

Returns `true`.
"""
issparse(S::SparseMatrixCSR) = true

"""
    nnz(S::SparseMatrixCSR)

Returns the number of stored (filled) elements in a sparse array.
"""
nnz(S::SparseMatrixCSR{Bi}) where Bi = Int(S.rowptr[size(S, 1) + 1] - Bi)

"""
    nonzeros(S::SparseMatrixCSR)

Return a vector (1-based) of the structural nonzero values in sparse array S.
This includes zeros that are explicitly stored in the sparse array.
The returned vector points directly to the internal nonzero storage of S,
and any modifications to the returned vector will mutate S as well.
"""
nonzeros(S::SparseMatrixCSR) = S.nzval

nzvalview(S::SparseMatrixCSR) = view(nonzeros(S), 1:nnz(S))

"""
    colvals(S::SparseMatrixCSR{Bi}) where {Bi}

Return a vector of the col indices of `S`. The stored values are indexes to arrays
with `Bi`-based indexing, but the `colvals(S)` array itself is a standard 1-based
Julia `Vector`.
Any modifications to the returned vector will mutate S as well.
Providing access to how the col indices are stored internally
can be useful in conjunction with iterating over structural
nonzero values. See also [`nonzeros`](@ref) and [`nzrange`](@ref).
"""
colvals(S::SparseMatrixCSR) = S.colval

"""
    nzrange(S::SparseMatrixCSR{Bi}, row::Integer) where {Bi}

Return the range of indices to the structural nonzero values of a
sparse matrix row. The returned range of indices is always 1-based even for `Bi != 1`.
"""
nzrange(S::SparseMatrixCSR{Bi}, row::Integer) where {Bi} = S.rowptr[row]+getoffset(S):S.rowptr[row+1]-Bi

"""
    findnz(S::SparseMatrixCSR{Bi,Tv,Ti})

Return a tuple `(I, J, V)` where `I` and `J` are the row and column 1-based indices
of the stored ("structurally non-zero") values in sparse matrix A,
and V is a vector of the values. The returned vectors are newly allocated
and are unrelated to the internal storage of matrix `S`.
"""
function findnz(S::SparseMatrixCSR{Bi,Tv,Ti}) where {Bi,Tv,Ti}
  numnz = nnz(S)
  I = Vector{Ti}(undef, numnz)
  J = Vector{Ti}(undef, numnz)
  V = Vector{Tv}(undef, numnz)
  count = 1
  o = getoffset(S)
  @inbounds for row in 1:size(S,1)
    @inbounds for k in nzrange(S,row)
      I[count] = row
      J[count] = S.colval[k]+o
      V[count] = S.nzval[k]
      count += 1
    end
  end
  return (I, J, V)
end

"""
    count(pred, S::SparseMatrixCSR)
    count(S::SparseMatrixCSR)

Count the number of elements in `nonzeros(S)` for which predicate `pred` returns `true`.
If  `pred` not given, it counts the number of `true` values.
"""
count(pred, S::SparseMatrixCSR) = count(pred, nzvalview(S))
count(S::SparseMatrixCSR) = count(i->true, nzvalview(S))

function mul!(y::AbstractVector,A::SparseMatrixCSR,v::AbstractVector, α::Number, β::Number)
  A.n == size(v, 1) || throw(DimensionMismatch())
  A.m == size(y, 1) || throw(DimensionMismatch())
  if β != 1
    β != 0 ? rmul!(y, β) : fill!(y, zero(eltype(y)))
  end
  o = getoffset(A)
  for row = 1:size(y, 1)
    @inbounds for nz in nzrange(A,row)
      col = A.colval[nz]+o
      y[row] += A.nzval[nz]*v[col]*α
    end
  end
  return y
end

function mul!(y::AbstractVector,A::SparseMatrixCSR,v::AbstractVector)
  A.n == size(v, 1) || throw(DimensionMismatch())
  A.m == size(y, 1) || throw(DimensionMismatch())
  fill!(y, zero(eltype(y)))
  o = getoffset(A)
  for row = 1:size(y, 1)
    @inbounds for nz in nzrange(A,row)
      col = A.colval[nz]+o
      y[row] += A.nzval[nz]*v[col]
    end
  end
  return y
end

*(A::SparseMatrixCSR, v::Vector) = (y = similar(v,size(A,1));mul!(y,A,v))

function show(io::IO, ::MIME"text/plain", S::SparseMatrixCSR)
  xnnz = nnz(S)
  print(io, S.m, "×", S.n, " ", typeof(S), " with ", xnnz, " stored ",
        xnnz == 1 ? "entry" : "entries")
  if xnnz != 0
    print(io, ":")
    show(IOContext(io, :typeinfo => eltype(S)), S)
  end
end
show(io::IO, S::SparseMatrixCSR) = show(convert(IOContext, io), S::SparseMatrixCSR)

function show(io::IOContext, S::SparseMatrixCSR{Bi}) where{Bi}
  nnz(S) == 0 && return show(io, MIME("text/plain"), S)
  ioc = IOContext(io, :compact => true)

  function _format_line(r, col, padr, padc, o=0)
    print(ioc, "\n  [", rpad(col+o, padr), ", ", lpad(S.colval[r]+o, padc), "]  =  ")
    if isassigned(S.nzval, Int(r))
      show(ioc, S.nzval[r])
    else
      print(ioc, Base.undef_ref_str)
    end
  end

  function _get_cols(from, to)
    idx = eltype(S.rowptr)[]
    c = searchsortedlast(S.rowptr, from)
    for i = from:to
      while i == S.rowptr[c+1]
        c +=1
      end
      push!(idx, c)
    end
    idx
  end

  rows = displaysize(io)[1] - 4 # -4 from [Prompt, header, newline after elements, new prompt]
  if !get(io, :limit, false) || rows >= nnz(S) # Will the whole matrix fit when printed?
    cols = _get_cols(Bi, nnz(S)-getoffset(S))
    padr, padc = ndigits.((maximum(S.colval[1:nnz(S)]), cols[end]))
    _format_line.(1:nnz(S), cols.-getoffset(S), padr, padc,getoffset(S))
  else
    if rows <= 2
      print(io, "\n  \u22ee")
      return
    end
    s1, e1 = 1, div(rows - 1, 2) # -1 accounts for \vdots
    s2, e2 = nnz(S) - (rows - 1 - e1) + 1, nnz(S)
    cols1, cols2 = _get_cols(s1-getoffset(S), e1-getoffset(S)), _get_cols(s2-getoffset(S), e2-getoffset(S))
    padr = ndigits(max(maximum(S.colval[s1:e1]), maximum(S.colval[s2:e2])))
    padc = ndigits(cols2[end])
    _format_line.(s1:e1, cols1, padr, padc)
    print(io, "\n  \u22ee")
    _format_line.(s2:e2, cols2, padr, padc)
  end
  return
end
