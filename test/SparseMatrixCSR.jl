@testset "SparseMatrixCSR" begin
    maxnz=5
    maxrows=5
    maxcols=5

    for T in (Int32,Int64,Float32,Float64)
        for Bi in (0,1)
            I = Vector{Int}()
            J = Vector{Int}()
            V = Vector{T}()
            for (ik, jk, vk) in zip(rand(1:maxrows, maxnz), rand(1:maxcols, maxnz), rand(1:T(maxnz), maxnz))
                push_coo!(SparseMatrixCSR,I,J,V,ik,jk,vk)
            end
            finalize_coo!(SparseMatrixCSR,I,J,V,maxcols,maxrows)
            CSC = sparse(I, J, V, maxcols,maxrows)
            CSR = sparsecsr(SparseMatrixCSR{Bi},I, J, V,maxcols,maxrows)

            @test CSC == CSR

            @test nnz(CSC) == count(!iszero, CSC) == nnz(CSR) == count(!iszero, CSR)


            TCSC = sparse(J, I, V, maxrows, maxcols)
            TCSR = sparsecsr(SparseMatrixCSR{Bi}, J, I, V, maxrows, maxcols)

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
    
end

