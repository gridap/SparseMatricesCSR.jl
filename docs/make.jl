using Documenter, SparseMatricesCSR

makedocs(;
    modules=[SparseMatricesCSR],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/gridap/SparseMatricesCSR.jl/blob/{commit}{path}#L{line}",
    sitename="SparseMatricesCSR.jl",
    authors=["Víctor Sande <vsande@cimne.upc.edu>","Francesc Verdugo <fverdugo@cimne.upc.edu>"],
    assets=String[],
)

deploydocs(;
    repo="github.com/gridap/SparseMatricesCSR.jl",
)
