using Documenter, SparseMatrix

makedocs(;
    modules=[SparseMatrix],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/gridap/SparseMatrix.jl/blob/{commit}{path}#L{line}",
    sitename="SparseMatrix.jl",
    authors="Víctor Sande <vsande@cimne.upc.edu>",
    assets=String[],
)

deploydocs(;
    repo="github.com/gridap/SparseMatrix.jl",
)
