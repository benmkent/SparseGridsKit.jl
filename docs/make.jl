using Documenter, SparseGridsKit

makedocs(
    sitename="SparseGridsKit",
    pages = [
    "SparseGridsKit" => "index.md",
    "Examples" => ["Multi-Index Sets"=>"multiindexsets.md",
                    "Sparse Grid Interpolation"=>"sparsegridinterpolation.md",
                    "Sparse Grid Integration"=>"sparsegridintegration.md"
                    ],
    "Reference" => "reference.md"
    ]
)

deploydocs(
    repo = "github.com/benmkent/SparseGridsKit.jl.git",
)
