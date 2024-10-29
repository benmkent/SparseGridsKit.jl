using Test, Documenter, SparseGridsKit

@testset "SparseGridsKit.jl tests" begin

    @testset "Sparse Grid tests" begin
        include("sparsegrids_test.jl")
    end

    @testset "doctest" begin
        doctest(SparseGridsKit)
    end
end