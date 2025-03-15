module SymSparseMatrixCSRTests

using Test
using SparseMatricesCSR
using SparseArrays
using LinearAlgebra

function test_csr(Bi,Tv,Ti)
  maxrows=5
  maxcols=6

  I = Ti[1,2,1,3,2,5,2,5]
  J = Ti[1,1,2,2,3,2,5,5]
  V = Tv[2,3,3,4,4,7,7,8]

  I_up = Ti[1,1,2,2,5]
  J_up = Ti[1,2,3,5,5]
  V_up = Tv[2,3,4,7,8]

  CSC = sparse(I,J,V)
  if Bi == 1
    CSR = symsparsecsr(I_up,J_up,V_up)
    @test CSR == CSC
    CSR = symsparsecsr(copy(I),copy(J),copy(V);symmetrize=true)
    @test CSR == CSC
    CSR = copy(symsparsecsr(copy(I),copy(J),copy(V);symmetrize=true))
    @test CSR == CSC
  end

  CSR = symsparsecsr(Val(Bi),I_up,J_up,V_up)
  @test CSR == CSC
  CSR = symsparsecsr(Val(Bi),copy(I),copy(J),copy(V);symmetrize=true)
  @test CSR == CSC
  @test eltype(CSR) == Tv
  @test isa(CSR,SymSparseMatrixCSR{Bi,Tv,Ti})

  for i=1:size(CSR,1)
    for j=1:size(CSR,2)
      if (i,j) in zip(I,J)
         CSR[i,j] = eltype(V)(i+j)
         @test CSR[i,j] ≈ eltype(V)(i+j)
      else
         try
          CSR[i,j] = eltype(V)(i+j)
         catch e
          @test isa(e,ArgumentError)
         end
      end
    end
  end

  CSC = sparse(I,J,V,maxrows,maxcols)
  if Bi == 1
    CSR = symsparsecsr(I_up,J_up,V_up,maxrows,maxcols)
    @test CSR == CSC
    CSR = symsparsecsr(copy(I),copy(J),copy(V),maxrows,maxcols;symmetrize=true)
    @test CSR == CSC
  end

  CSR = symsparsecsr(Val(Bi),I_up,J_up,V_up,maxrows,maxcols)
  @test CSR == CSC
  CSR = symsparsecsr(Val(Bi),copy(I),copy(J),copy(V),maxrows,maxcols;symmetrize=true)
  @test CSR == CSC

  CSR2 = convert(typeof(CSR),CSR)
  @test CSR2 === CSR
  CSR2 = convert(SymSparseMatrixCSR{0,Float64,Int32},CSR)
  @test CSR2 == CSR

  @test size(CSR) == size(CSC)
  @test eltype(CSR) == Tv
  @test isa(CSR,SymSparseMatrixCSR{Bi,Tv,Ti})
  @test issparse(CSR)
  @test getBi(CSR) == Bi
  @test getoffset(CSR) == 1-Bi
  @test nnz(CSR) == nnz(CSR.uppertrian)
  @test length(nonzeros(CSR)) == nnz(CSR)
  @test nonzeros(CSR) === CSR.uppertrian.nzval
  @test colvals(CSR) === CSR.uppertrian.colval
  i,j,v = findnz(CSR)
  csr = symsparsecsr(Val(Bi),i,j,v,maxrows,maxcols)
  @test csr == CSR
  @test count(v->v>0,CSR) == count(v->v>0,nonzeros(CSR))

  x = rand(Tv,maxcols)
  y = Vector{Tv}(undef,maxrows)
  z = Vector{Tv}(undef,maxrows)
  mul!(y,CSR,x)
  mul!(z,CSC,x)
  @test y ≈ z

  mul!(y,CSR,x,1,2)
  mul!(z,CSC,x,1,2)
  @test y ≈ z

  @test CSR*x ≈ CSC*x

  out = LinearAlgebra.fillstored!(CSR,3.33)
  @test out === CSR
  LinearAlgebra.fillstored!(CSC,3.33)
  mul!(y,CSR,x)
  mul!(z,CSC,x)
  @test y ≈ z

  _CSR = copy(CSR)
  out = LinearAlgebra.rmul!(CSR,-1)
  @test out === CSR
  @test _CSR ≈ -1*CSR

end

for Bi in (0,1)
  for Tv in (Float32,Float64)
    for Ti in (Int32,Int64)
      test_csr(Bi,Tv,Ti)
    end
  end
end

end # module
