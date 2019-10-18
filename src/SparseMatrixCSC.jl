
"""
    function push_coo!(::Type{SparseMatrixCSC},I,J,V,ik,jk,vk) 

Inserts entries in COO vectors for further building a SparseMatrixCSC.
"""
function push_coo!(::Type{SparseMatrixCSC},I,J,V,ik,jk,vk) 
    (push!(I, ik), push!(J, jk), push!(V, vk))
end
