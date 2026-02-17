using Documenter, DocumenterInterLinks
using SparseMatricesCSR

links = InterLinks(
    "SparseArrays" => "https://docs.julialang.org/en/v1/"
)

makedocs(;
    modules=[SparseMatricesCSR],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/gridap/SparseMatricesCSR.jl/blob/{commit}{path}#L{line}",
    sitename="SparseMatricesCSR.jl",
    authors="Víctor Sande <vsande@cimne.upc.edu> and Francesc Verdugo <fverdugo@cimne.upc.edu>",
    plugins=[links],
)

deploydocs(;
    repo="github.com/gridap/SparseMatricesCSR.jl",
)
