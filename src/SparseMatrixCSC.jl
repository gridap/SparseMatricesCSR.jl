"""
    function push_coo!(::Type{SparseMatrixCSC},I,J,V,ik,jk,vk) 

Inserts entries in COO vectors for further building a SparseMatrixCSC.
"""
function push_coo!(::Type{SparseMatrixCSC},I::Vector{Integer},J::Vector{Integer},V::Vector{Number},ik::Integer,jk::Integer,vk::Number)
    (push!(I, ik), push!(J, jk), push!(V, vk))
end

"""
    function push_coo!(I,J,V,ik,jk,vk) 

Inserts entries in COO vectors for further building a SparseMatrixCSC.
"""
function push_coo!(I::Vector{Integer},J::Vector{Integer},V::Vector{Number},ik::Integer,jk::Integer,vk::Number)
    push_coo!(SparseMatrixCSC,J,V,ik,jk,vk)
end

"""
    function finalize_coo!(::Type{SparseMatrixCSC},I,J,V,m,n) 

Check and insert diagonal entries in COO vectors if needed.
"""
function finalize_coo!(::Type{SparseMatrixCSC},I::Vector{Integer},J::Vector{Integer},V::Vector{Number},m::Integer,n::Integer)
end

"""
    function finalize_coo!(I,J,V,m,n) 

Check and insert diagonal entries in COO vectors if needed.
"""
function finalize_coo!(I::Vector{Integer},J::Vector{Integer},V::Vector{Number},m::Integer,n::Integer)
    finalize_coo!(SparseMatrixCSC,I,J,V,m,n)
end
