@testset "SparseMatrixCSC" begin
    maxnz=5
    maxrows=5
    maxcols=5

    for T in (Int32,Int64,Float32,Float64)
        I = Vector{Int}()
        J = Vector{Int}()
        V = Vector{T}()
        for (ik, jk, vk) in zip(rand(1:maxrows, maxnz), rand(1:maxcols, maxnz), rand(1:T(maxnz), maxnz))
            push_coo!(SparseMatrixCSC,I,J,V,ik,jk,vk)
        end
        finalize_coo!(SparseMatrixCSC,I,J,V,maxcols,maxrows)
        CSC = sparse(I, J, V, maxcols,maxrows)

        for col = 1:maxcols
            for j in nzrange(CSC, col)
                @test nonzeros(CSC)[j] == CSC[rowvals(CSC)[j], col]
            end
        end

        @test size(CSC) == (maxrows,maxcols)

        @test nnz(CSC) == count(i->(i!=0), CSC) <= maxnz

    end
    
end
