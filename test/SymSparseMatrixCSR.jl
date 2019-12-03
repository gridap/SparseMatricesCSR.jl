@testset "SymSparseMatrixCSR" begin
    maxnz=5
    maxrows=5
    maxcols=5

    for T in (Int32,Int64,Float32,Float64)
        for Bi in (0,1)
            I = Vector{Int}()
            J = Vector{Int}()
            V = Vector{T}()
            for (ik, jk, vk) in zip(rand(1:maxrows, maxnz), rand(1:maxcols, maxnz), rand(1:T(maxnz), maxnz) )
                push_coo!(SymSparseMatrixCSR,I,J,V,ik,jk,vk)
            end
            finalize_coo!(SymSparseMatrixCSR,I,J,V,maxrows, maxcols)
            SYMCSC = Symmetric(sparse(I, J, V, maxrows, maxcols),:U)
            SYMCSR = symsparsecsr(SymSparseMatrixCSR{Bi},I, J, V, maxrows, maxcols)

            @test size(SYMCSC)==size(SYMCSR)
            @test SYMCSC == SYMCSR

            @test convert(SymSparseMatrixCSR{Bi}, SYMCSR) === SYMCSR

            SYMCSRC = convert(SymSparseMatrixCSR{SYMCSR.uppertrian.offset}, SYMCSR)
            @test SYMCSRC == SYMCSR
            @test SYMCSRC !== SYMCSR

            @test hasrowmajororder(SYMCSR) == true
            @test hascolmajororder(SYMCSR) == false
            @test getptr(SYMCSR)           == SYMCSR.uppertrian.rowptr
            @test getindices(SYMCSR)       == colvals(SYMCSR)

            @test nnz(SYMCSC.data) == nnz(SYMCSR.uppertrian) <= nnz(SYMCSR) 
            @test count(!iszero, SYMCSC.data) == count(!iszero, SYMCSR.uppertrian)

            ICSC,JCSC,VCSC= findnz(SYMCSC.data)
            ICSR,JCSR,VCSR= findnz(SYMCSR)

            @test sort(ICSC)==sort(JCSR) && sort(JCSC)==sort(ICSR) && sort(VCSC)==sort(VCSR)

            v = rand(size(SYMCSC)[2])
            @test SYMCSC*v == SYMCSR*v
        end
    end
    
end

