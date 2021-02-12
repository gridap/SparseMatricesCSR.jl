module SparseMatrixCSRTests

using Test
using SparseMatricesCSR
using SparseArrays
using LinearAlgebra

function test_csr(Bi,Tv,Ti)
  maxnz=10
  maxrows=5
  maxcols=5

  I = rand(Ti(1):Ti(maxrows),maxnz)
  J = rand(Ti(1):Ti(maxcols),maxnz)
  V = rand(Tv,maxnz)

  #I = Ti[1,2,3]
  #J = Ti[1,1,3]
  #V = Tv[3,4,6]

  CSC = sparse(I,J,V)
  if Bi == 1
    CSR = sparsecsr(I,J,V)
    @test CSR == CSC
  end
  CSR = sparsecsr(Val(Bi),I,J,V)
  @test CSR == CSC
  CSR = sparsecsr(Val(Bi),I,J,V)
  @test CSR == CSC
  @test eltype(CSR) == Tv
  @test isa(CSR,SparseMatrixCSR{Bi,Tv,Ti})

  CSC = sparse(I,J,V,maxrows,maxcols)
  if Bi == 1
    CSR = sparsecsr(I,J,V,maxrows,maxcols)
    @test CSR == CSC
  end
  CSR = sparsecsr(Val(Bi),I,J,V,maxrows,maxcols)
  @test CSR == CSC
  CSR = sparsecsr(Val(Bi),I,J,V,maxrows,maxcols)
  @test CSR == CSC

  @test size(CSR) == size(CSC)
  @test eltype(CSR) == Tv
  @test isa(CSR,SparseMatrixCSR{Bi,Tv,Ti})
  @test issparse(CSR)
  @test getBi(CSR) == Bi
  @test getoffset(CSR) == 1-Bi
  @test nnz(CSR) == nnz(CSC)
  @test length(nonzeros(CSR)) == nnz(CSR)
  @test nonzeros(CSR) === CSR.nzval
  @test colvals(CSR) === CSR.colval
  i,j,v = findnz(CSR)
  csr = sparsecsr(Val(Bi),i,j,v,maxrows,maxcols)
  @test csr == CSR
  @test count(v->v>0,CSR) == count(v->v>0,nonzeros(CSR))

  x = rand(Tv,maxcols)
  y = Vector{Tv}(undef,maxcols)
  z = Vector{Tv}(undef,maxcols)
  mul!(y,CSR,x)
  mul!(z,CSC,x)
  @test y ≈ z

  mul!(y,CSR,x,1,2)
  mul!(z,CSC,x,1,2)
  @test y ≈ z

  @test CSR*x ≈ CSC*x

end

for Bi in (0,1)
  for Tv in (Float32,Float64)
    for Ti in (Int32,Int64)
      test_csr(Bi,Tv,Ti)
    end
  end
end

end # module


#@testset "SparseMatrixCSR" begin
#    maxnz=5
#    maxrows=5
#    maxcols=5
#    int_types=(Int32,Int64)
#    float_types=(Float32,Float64)
#    Bi_types=(0,1)
#
#    for Ti in int_types
#        for Tv in float_types
#            for Bi in Bi_types
#                I = Vector{Ti}()
#                J = Vector{Ti}()
#                V = Vector{Tv}()
#                for (ik, jk, vk) in zip(rand(1:maxrows, maxnz), rand(1:maxcols, maxnz), rand(1:Tv(maxnz), maxnz))
#                    push_coo!(SparseMatrixCSR{Bi,Tv,Ti},I,J,V,ik,jk,vk)
#                end
#                finalize_coo!(SparseMatrixCSR{Bi,Tv,Ti},I,J,V,maxrows,maxcols)
#                CSC = sparse(I, J, V, maxrows,maxcols)
#                CSR = sparsecsr(SparseMatrixCSR{Bi,Tv,Ti},I, J, V,maxrows,maxcols)
#
#                show(devnull, CSR);
#
#                @test CSC == CSR
#                @test nnz(CSC) == count(!iszero, CSC) == nnz(CSR) == count(!iszero, CSR)
#
#                @test hasrowmajororder(CSR) == true
#                @test hascolmajororder(CSR) == false
#                @test getptr(CSR)           == CSR.rowptr
#                @test getindices(CSR)       == colvals(CSR)
#
#                TCSC = sparse(J, I, V, maxrows, maxcols)
#                TCSR = sparsecsr(SparseMatrixCSR{Bi}, J, I, V, maxrows, maxcols)
#
#                @test size(CSC)==size(CSR)==reverse(size(TCSC))==reverse(size(TCSC))
#
#                @test [nzrange(CSC,col) for col in 1:size(CSC,2)] == [nzrange(TCSR,row) for row in 1:size(TCSR,1)]
#                @test [nzrange(CSR,row) for row in 1:size(CSR,1)] == [nzrange(TCSC,col) for col in 1:size(TCSC,2)]
#
#                @test nonzeros(CSC) == nonzeros(TCSR) && nonzeros(CSR) == nonzeros(TCSC) 
#
#                ICSC,JCSC,VCSC= findnz(CSC)
#                ICSR,JCSR,VCSR= findnz(CSR)
#
#                @test sort(ICSC)==sort(JCSR) && sort(JCSC)==sort(ICSR) && sort(VCSC)==sort(VCSR)
#
#                v = rand(size(CSC)[2])
#                @test CSC*v == CSR*v
#
#                for cBi in Bi_types
#                    if Bi == cBi
#                        @test convert(SparseMatrixCSR{Bi}, CSR) === CSR
#                    else
#                        CSRC = convert(SparseMatrixCSR{cBi}, CSR)
#                        @test CSRC == CSR
#                        @test CSRC !== CSR
#                    end
#                end
#
#                for cTi in int_types
#                    for cTv in float_types
#                        if (Ti,Tv) == (cTi,cTv)
#                            @test convert(SparseMatrixCSR{Bi,Tv,Ti}, CSR) === CSR
#                        else
#                            CSRC = convert(SparseMatrixCSR{Bi,cTv,cTi}, CSR)
#                            @test CSRC == CSR
#                            @test CSRC !== CSR
#                        end
#                    end
#                end
#
#            end
#        end
#    end    
#end

