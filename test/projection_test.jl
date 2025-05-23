using Test
using SparseGridsKit

@testset "Projection" begin
    # Test case 1: Simple projection
    sg = create_sparsegrid(create_smolyak_miset(1, 5))
    n=10
    x = [[xi] for xi = 0:0.1:1]
    f(x) = x[1]^3 + x[1]^2+3 
    y = [[f(xi)] for xi in x]

    sga = project_onto_sparsegrid(sg, x, y)
    @test all(sga(x) ≈ [f(x)] for x in 0:0.05:1)


    # Test case 2: Edge case with empty input
    x_empty = []
    y_empty = []
    result_empty = project_onto_sparsegrid(sg, x_empty, y_empty)
    @test length(result_empty[2]) == 0
end