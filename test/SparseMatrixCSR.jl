module SparseMatrixCSRTests

using Test
using SparseMatricesCSR
using SparseArrays
using LinearAlgebra

function test_csr(Bi,Tv,Ti)
  maxnz=10
  maxrows=5
  maxcols=6

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
    @test copy(CSR) == CSC
  end
  CSR = sparsecsr(Val(Bi),I,J,V)
  show(IOContext(stdout, :limit=>true, :displaysize=>(10,10)), CSR)
  show(IOContext(stdout, :limit=>false), CSR)
  @test CSR == CSC
  @test copy(CSR) == CSC
  @test eltype(CSR) == Tv
  @test isa(CSR,SparseMatrixCSR{Bi,Tv,Ti})

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

  if Ti == Int64
      dense = rand(Tv, maxrows,maxcols)
      CSR = sparsecsr(Val(Bi), dense)
      @test isa(CSR,SparseMatrixCSR{Bi,Tv,Ti})
      @test CSR == dense
      if Bi == 1
          CSR = sparsecsr(dense)
          @test CSR == dense
          @test isa(CSR,SparseMatrixCSR{Bi,Tv,Ti})
      end
  end

  CSC = sparse(I,J,V,maxrows,maxcols)
  if Bi == 1
    CSR = sparsecsr(I,J,V,maxrows,maxcols)
    @test CSR == CSC
    @test copy(CSR) == CSC
  end
  CSR = sparsecsr(Val(Bi),I,J,V,maxrows,maxcols)
  @test CSR == CSC
  @test copy(CSR) == CSC
  CSR = sparsecsr(Val(Bi),I,J,V,maxrows,maxcols)
  @test CSR == CSC
  @test copy(CSR) == CSC

  CSR2 = convert(typeof(CSR),CSR)
  @test CSR2 === CSR
  CSR2 = convert(SparseMatrixCSR{0,Float64,Int32},CSR)
  @test CSR2 == CSR
  @test copy(CSR2) == CSR
  CSR2 = convert(SparseMatrixCSR{0,Float64,Int32},CSC)
  @test CSR2 == CSC
  @test copy(CSR2) == CSC
  CSR2 = convert(SparseMatrixCSR{0,Float64,Int32},collect(CSC))
  @test CSR2 == CSC
  @test copy(CSR2) == CSC

  @test size(CSR) == size(CSC)
  @test eltype(CSR) == Tv
  @test isa(CSR,SparseMatrixCSR{Bi,Tv,Ti})
  @test issparse(CSR)
  @test getBi(CSR) == Bi
  @test getoffset(CSR) == 1-Bi
  @test nnz(CSR) == nnz(CSC)
  @test length(SparseArrays.nzvalview(CSR)) == nnz(CSR)
  @test nonzeros(CSR) === CSR.nzval
  @test colvals(CSR) === CSR.colval
  i,j,v = findnz(CSR)
  csr = sparsecsr(Val(Bi),i,j,v,maxrows,maxcols)
  @test csr == CSR
  @test count(v->v>0,CSR) == count(v->v>0,nonzeros(CSR))

  x = rand(Tv,maxcols)
  y = Vector{Tv}(undef,maxrows)
  z = Vector{Tv}(undef,maxrows)
  mul!(y,CSR,x)
  mul!(z,CSC,x)
  @test y ≈ z
  @test CSR*x ≈ CSC*x

  mul!(y,CSR,x,1,2)
  mul!(z,CSC,x,1,2)
  @test y ≈ z

  out = LinearAlgebra.fillstored!(CSR,3.33)
  @test out === CSR
  LinearAlgebra.fillstored!(CSC,3.33)
  mul!(y,CSR,x)
  mul!(z,CSC,x)
  @test y ≈ z

  # test constructors
  @test CSR == SparseMatrixCSR(CSC)
  @test CSR == SparseMatrixCSR(Matrix(CSC))

  _CSR = copy(CSR)
  out = LinearAlgebra.rmul!(CSR,-1)
  @test out === CSR
  @test _CSR ≈ -1*CSR

  # Test overlong buffers
  let A = sparsecsr([1, 2, 3], [1, 2, 3], [4., 5, 6])
    push!(A.nzval, 4.0)
    push!(A.rowptr, -1)
    @test nnz(A) == length(SparseArrays.nzvalview(A)) == 3
    @test SparseArrays.nzvalview(A) == [4., 5, 6]
  end
end

function test_lu(Bi,I,J,V)
  CSR=sparsecsr(Val(Bi),I,J,V)
  CSC=sparse(I,J,V)
  x=rand(3)
  @test norm(CSR\x-CSC\x) < 1.0e-14
  fact=lu(CSR)
  lu!(fact,CSR)
  y=similar(x)
  ldiv!(y,fact,x)
  @test norm(y-CSC\x) < 1.0e-14
end


for Bi in (0,1)
  for Tv in (Float32,Float64)
    for Ti in (Int32,Int64)
      test_csr(Bi,Tv,Ti)
    end
  end
end

if Base.USE_GPL_LIBS  # `lu!` requires `SuiteSparse.UMFPACK`
  I = [1,1,2,2,2,3,3]
  J = [1,2,1,2,3,2,3]
  V = [4.0,1.0,-1.0,4.0,1.0,-1.0,4.0]
  test_lu(0,I,J,V)
  test_lu(1,I,J,V)
else
  @warn "Tests run without GPL libraries."
end

end # module
