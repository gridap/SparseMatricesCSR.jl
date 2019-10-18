@testset "SymSparseMatrixCSR" begin
    maxnz=5
    maxrows=5
    maxcols=5

    for T in (Int32,Int64,Float32,Float64)
        I = vcat(rand(1:maxrows, maxnz),1:maxrows)
        J = vcat(rand(1:maxcols, maxnz),1:maxcols)
        V = vcat(rand(1:T(maxnz), maxnz), zeros(T,maxnz))
        indices = [ i for i in 1:length(I) if !(I[i]<J[i])] # Lower triangle indices
        SYMCSC = Symmetric(sparse(I[indices], J[indices], V[indices], maxrows, maxcols),:L)
        SYMCSR = symsparsecsr(I[indices], J[indices], V[indices], maxrows, maxcols)

        @test size(SYMCSC)==size(SYMCSR)
        @test SYMCSC == SYMCSR

        @test nnz(SYMCSC.data) == nnz(SYMCSR.lowertrian) <= nnz(SYMCSR) 
        @test count(i->(i!=0), SYMCSC.data) == count(i->(i!=0), SYMCSR.lowertrian)

        ICSC,JCSC,VCSC= findnz(SYMCSC.data)
        ICSR,JCSR,VCSR= findnz(SYMCSR)

        @test sort(ICSC)==sort(JCSR) && sort(JCSC)==sort(ICSR) && sort(VCSC)==sort(VCSR)

        v = rand(size(SYMCSC)[2])
        @test SYMCSC*v == SYMCSR*v
    end
    
end

