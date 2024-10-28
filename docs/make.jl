using Documenter, SparseGridsKit

makedocs(
    sitename="SparseGridsKit",
    pages = [
    "SparseGridsKit" => "index.md",
    "Examples" => ["Multi-Index Sets"=>"multiindexsets.md"]
    ]
)