
struct SparseMatrixCSR{T,Ti<:Integer} <: AbstractSparseMatrix{T,Ti}
    transpose :: SparseMatrixCSC{T,Ti}
end


show(io::IO, A::SparseMatrixCSR) = show(io, A.transpose)
size(A::SparseMatrixCSR) = (A.transpose.n, A.transpose.m)
getindex(A::SparseMatrixCSR, x::Integer, y::Integer) = getindex(A.transpose,y,x)

nnz(S::SparseMatrixCSR)         = nnz(S.transpose)
count(pred, S::SparseMatrixCSR) = count(pred, S.transpose)

nonzeros(S::SparseMatrixCSR) = S.transpose.nzval

nzrange(S::SparseMatrixCSR, row::Integer) = nzrange(S.transpose, row)

findnz(S::SparseMatrixCSR) = findnz(S.transpose)

sparsecsr(I,J,kwargs...) = SparseMatrixCSR(sparse(J,I,kwargs...))


