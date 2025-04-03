using SparseGridsKit, Documenter, DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))
makedocs(
    sitename="SparseGridsKit.jl",
    pages = [
    "SparseGridsKit.jl" => "index.md",
    "SparseGridsKit.jl by Examples" => ["Knots"=>"knots.md",
                    "Multi-Index Sets"=>"misets.md",
                    "Sparse Grids"=>"sparsegrids.md",
                    "Sparse Grid Interpolation"=>"sparsegridinterpolation.md",
                    "Sparse Grid Integration"=>"sparsegridintegration.md",
                    "Adaptive Sparse Grids"=>"adaptivesparsegrids.md",
                    "Spectral Sparse Grids"=>"spectralsparsegrids.md",
                    ],
    "Plotting" => "plots.md",
    "Extensions" => ["Interpolation of PDE solutions"=>"functioninterpolation.md",
    "UM-Bridge" => "umbridge.md"],
    "SparseGridsKit.jl Public API Reference"=>"reference.md",
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