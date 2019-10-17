@testset "SymSparseMatrixCSR" begin
    maxnz=5
    maxrows=5
    maxcols=5

    for T in (Int32,Int64,Float32,Float64)
        I = vcat(rand(1:maxrows, maxnz),1:maxrows)
        J = vcat(rand(1:maxcols, maxnz),1:maxcols)
        V = vcat(rand(1:T(maxnz), maxnz), zeros(T,maxnz))
        SYMCSC = Symmetric(sparse(I, J, V, maxrows, maxcols),:L)
        SYMCSR = symsparsecsr(I, J, V, maxrows, maxcols)

        @test size(SYMCSC)==size(SYMCSR)

        @test SYMCSC == SYMCSR

        v = rand(size(SYMCSC)[2])
        @test SYMCSC*v == SYMCSR*v
    end
    
end

