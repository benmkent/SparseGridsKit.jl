using Test

@testset "SparseGridsKit.jl tests" begin

    @testset "Sparse Grid tests" begin
        include("sparsegrids_tests.jl")
    end

end