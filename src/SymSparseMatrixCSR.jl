

struct SymSparseMatrixCSR{T,Ti<:Integer} <: AbstractSparseMatrix{T,Ti}
    lowertrian :: SparseMatrixCSC{T,Ti}
end


show(io::IO, A::SymSparseMatrixCSR) = show(io, A.lowertrian)
size(A::SymSparseMatrixCSR) = size(A.lowertrian)
getindex(A::SymSparseMatrixCSR, x::Integer, y::Integer) = getindex(A.lowertrian,min(x,y),max(x,y))

nnz(S::SymSparseMatrixCSR) = 2*(nnz(S.lowertrian)-S.lowertrian.m)+S.lowertrian.m
count(pred, S::SymSparseMatrixCSR) = count(pred, S.lowertrian)
nonzeros(S::SymSparseMatrixCSR) = S.lowertrian.nzval
nzrange(S::SymSparseMatrixCSR, row::Integer) = nzrange(S.lowertrian, row)
findnz(S::SymSparseMatrixCSR) = findnz(S.lowertrian)
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

