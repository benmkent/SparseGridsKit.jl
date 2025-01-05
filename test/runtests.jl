using Test, Documenter, SparseGridsKit, LinearAlgebra

@testset "SparseGridsKit.jl tests" begin
    include("multiindexsets_test.jl")
    include("integration_test.jl")
    include("interpolation_test.jl")
    include("sparsegrids_test.jl")
    include("mixedpoints_test.jl")
    include("adaptive_test.jl")
end

@testset "doctest" begin
    doctest(SparseGridsKit)
end