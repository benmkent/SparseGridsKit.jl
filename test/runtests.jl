using Test, SparseGridsKit, Documenter, LinearAlgebra

@testset "SparseGridsKit.jl tests" begin
    include("knots_test.jl")
    include("multiindexsets_test.jl")
    include("sparsegrids_test.jl")
    include("interpolation_test.jl")
    include("adaptive_test.jl")
    include("multifidelity_tests.jl")
    include("integration_test.jl")
    include("spectral_test.jl")
    include("sgmk_test.jl")
    include("aqua_test.jl")
    include("plots_test.jl")
    include("differentiation_test.jl")
end