using Documenter, SparseGridsKit

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))
makedocs(
    sitename="SparseGridsKit",
    pages = [
    "SparseGridsKit" => "index.md",
    "Examples" => ["Multi-Index Sets"=>"multiindexsets.md",
                    "Sparse Grids"=>"sparsegrids.md",
                    "Sparse Grid Interpolation"=>"sparsegridinterpolation.md",
                    "Sparse Grid Integration"=>"sparsegridintegration.md"
                    ],
    "Reference" => "reference.md",
    "References" => "bib.md"
    ];
    format = Documenter.HTML(
        # ...
        assets=String["assets/citations.css"],
    ),
    plugins=[bib]
)

deploydocs(
    repo = "github.com/benmkent/SparseGridsKit.jl.git",
)
