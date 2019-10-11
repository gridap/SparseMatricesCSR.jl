@testset "SparseMatrixCSR" begin
    maxnz=5
    maxrows=5
    maxcols=5

    for T in (Int32,Int64,Float32,Float64)
        I = rand(1:maxrows, maxnz)
        J = rand(1:maxcols, maxnz)
        V = rand(1:T(maxnz),   maxnz)
        CSC = sparse(I, J, V)
        CSR = sparsecsr(I, J, V)

        @test size(CSC)==size(CSR)==reverse(size(CSR.mat))

        @test CSC==CSR

        @test nnz(CSC) == count(i->(i!=0), CSC) == nnz(CSR) == count(i->(i!=0), CSR)

        ICSC,JCSC,VCSC= findnz(CSC)
        ICSR,JCSR,VCSR= findnz(CSR)

        @test sort(ICSC)==sort(JCSR) && sort(JCSC)==sort(ICSR) && sort(VCSC)==sort(VCSR)

        v = rand(size(CSC)[2])
        @test CSC*v == CSR*v
    end
    
end

