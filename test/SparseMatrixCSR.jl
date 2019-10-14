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



        @test CSC==CSR

        @test nnz(CSC) == count(i->(i!=0), CSC) == nnz(CSR) == count(i->(i!=0), CSR)

        TCSC = sparse(J, I, V)
        TCSR = sparsecsr(J, I, V)

        @test size(CSC)==size(CSR)==reverse(size(CSR.transpose))
        @test size(CSC)==size(CSR)==reverse(size(TCSC))==reverse(size(TCSC))

        @test [nzrange(CSC,col) for col in 1:size(CSC,2)] == [nzrange(TCSR,row) for row in 1:size(TCSR,1)]
        @test [nzrange(CSR,row) for row in 1:size(CSR,1)] == [nzrange(TCSC,col) for col in 1:size(TCSC,2)]

        @test nonzeros(CSC) == nonzeros(TCSR) && nonzeros(CSR) == nonzeros(TCSC) 

        ICSC,JCSC,VCSC= findnz(CSC)
        ICSR,JCSR,VCSR= findnz(CSR)

        @test sort(ICSC)==sort(JCSR) && sort(JCSC)==sort(ICSR) && sort(VCSC)==sort(VCSR)

        v = rand(size(CSC)[2])
        @test CSC*v == CSR*v
    end
    
end

