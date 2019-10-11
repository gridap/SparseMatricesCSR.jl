
mutable struct SparseMatrixCSR{T,Ti<:Integer} <: AbstractSparseMatrix{T,Ti}
    mat :: SparseMatrixCSC{T,Ti}
    function SparseMatrixCSR{T,Ti}(m::Integer, n::Integer, rowptr::Vector{Ti}, colval::Vector{Ti}, nzval::Vector{T}) where {T,Ti<:Integer}
        csr = new()
        csr.mat = SparseMatrixCSC(m, n, rowptr, colval, nzval)
        return csr
    end
end

function SparseMatrixCSR(m::Integer, n::Integer, rowptr::Vector, colval::Vector, nzval::Vector)
    Tv = eltype(nzval)
    Ti = promote_type(eltype(rowptr), eltype(colval))
    return SparseMatrixCSR{Tv,Ti}(m, n, rowptr, colval, nzval)
end

const SparseMatrixCSRView{Tv,Ti} =
    SubArray{Tv,2,SparseMatrixCSR{Tv,Ti},
        Tuple{Base.Slice{Base.OneTo{Int}},I}} where {I<:AbstractUnitRange}
const SparseMatrixCSRUnion{Tv,Ti} = Union{SparseMatrixCSR{Tv,Ti}, SparseMatrixCSRView{Tv,Ti}}

show(io::IO, A::SparseMatrixCSR) = show(io, A.mat)
size(A::SparseMatrixCSR) = (A.mat.n, A.mat.m)
getindex(A::SparseMatrixCSR, x::Integer, y::Integer) = getindex(A.mat,y,x)

nnz(S::SparseMatrixCSR)         = nnz(S.mat)
count(pred, S::SparseMatrixCSR) = count(pred, S.mat)

nonzeros(S::SparseMatrixCSR) = S.mat.nzval
nonzeros(S::SparseMatrixCSRView) = nonzeros(S.mat.parent)

nzrange(S::SparseMatrixCSR, row::Integer) = nzrange(S.mat, row)
nzrange(S::SparseMatrixCSRView, row::Integer) = nzrange(S.mat.parent, S.indices[2][row])

findnz(S::SparseMatrixCSR) = findnz(S.mat)

function sparsecsr(I,J,kwargs...) 
    CSC=sparse(J,I,kwargs...)
    return SparseMatrixCSR(CSC.m, CSC.n, CSC.colptr, CSC.rowval, CSC.nzval)
end


