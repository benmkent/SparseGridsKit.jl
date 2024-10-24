using Test
using SparseGridsKit
using LinearAlgebra

@testset "Multi-Index Set Tests" begin
    # Test create_smolyak_miset
    n, k = 4, 3
    mi_set = create_smolyak_miset(n, k)
    @test get_n_mi(mi_set) == 35
    @test get_mi(mi_set) == [   [1, 1, 1, 1],
                                [1, 1, 1, 2],
                                [1, 1, 1, 3],
                                [1, 1, 1, 4],
                                [1, 1, 2, 1],
                                [1, 1, 2, 2],
                                [1, 1, 2, 3],
                                [1, 1, 3, 1],
                                [1, 1, 3, 2],
                                [1, 1, 4, 1],
                                [1, 2, 1, 1],
                                [1, 2, 1, 2],
                                [1, 2, 1, 3],
                                [1, 2, 2, 1],
                                [1, 2, 2, 2],
                                [1, 2, 3, 1],
                                [1, 3, 1, 1],
                                [1, 3, 1, 2],
                                [1, 3, 2, 1],
                                [1, 4, 1, 1],
                                [2, 1, 1, 1],
                                [2, 1, 1, 2],
                                [2, 1, 1, 3],
                                [2, 1, 2, 1],
                                [2, 1, 2, 2],
                                [2, 1, 3, 1],
                                [2, 2, 1, 1],
                                [2, 2, 1, 2],
                                [2, 2, 2, 1],
                                [2, 3, 1, 1],
                                [3, 1, 1, 1],
                                [3, 1, 1, 2],
                                [3, 1, 2, 1],
                                [3, 2, 1, 1],
                                [4, 1, 1, 1]];

    # Test add_mi
    mi_set_new = MISet([[1,1,1,1]]) 
    combined_mi_set = add_mi(mi_set, mi_set_new)
    @test get_n_mi(combined_mi_set) == 35
    mi_set_new = MISet([[5,1,1,1]]) 
    combined_mi_set = add_mi(mi_set, mi_set_new)
    @test get_n_mi(combined_mi_set) == 36

    # Test get_margin and get_reduced_margin
    mi_set = MISet([[1,1,1,1]]) 
    margin_set = get_margin(mi_set)
    @test get_mi(margin_set) ==  [  [1, 1, 1, 2],
                                    [1, 1, 2, 1],
                                    [1, 2, 1, 1],
                                    [2, 1, 1, 1]];
    
    reduced_margin_set = get_reduced_margin(mi_set)
    @test get_mi(reduced_margin_set) ==  [  [1, 1, 1, 2],
                                            [1, 1, 2, 1],
                                            [1, 2, 1, 1],
                                            [2, 1, 1, 1]];
end

@testset "Sparse Grids Tests" begin
    # Test create_sparsegrid
    n,k =4,0
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    @test get_n_grid_points(sg) == 1
    k=1
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    @test get_n_grid_points(sg) == 9
    k=2
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    @test get_n_grid_points(sg) == 41

    # Test get_grid_points and get_n_grid_points
    n,k =4,0
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    grid_points = get_grid_points(sg)
    @test length(grid_points) == 1
    @test grid_points[1] ≈ [0.0,0.0,0.0,0.0]
    n,k =1,3
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    grid_points = get_grid_points(sg)
    @test grid_points ≈ [   [0.0],
                            [1.0],
                            [-1.0],
                            [0.7071067811865475],
                            [-0.7071067811865475],
                            [0.9238795325112867],
                            [0.3826834323650898],
                            [-0.3826834323650898],
                            [-0.9238795325112867]]

    # Test map_from_to
    sg_from = create_sparsegrid(create_smolyak_miset(1, 8))
    sg_to = create_sparsegrid(create_smolyak_miset(1, 9))
    map_vector = map_from_to(sg_from, sg_to)
    @test map_vector == 1:get_n_grid_points(sg_from)
end

@testset "Sparse Grid Interpolation Test" begin
    # Test interpolate_on_sparsegrid
    f(x) = @. 3.0*x[1]^2 + 2.0*x[1] + 1.0

    n,k = 2,0
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    f_on_grid = [[f(x), cos(x[2])] for x in get_grid_points(sg)]
    target_points = [[x,x] for x in -1.0:0.1:1.0]
    interpolation_result = interpolate_on_sparsegrid(sg, f_on_grid, target_points)
    @test all(v ≈ [1.0, 1.0] for v in interpolation_result)

    n,k = 2,1
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    f_on_grid = [[f(x), cos(x[2])] for x in get_grid_points(sg)]
    target_points = [[x,x] for x in -1.0:0.1:1.0]
    interpolation_result = interpolate_on_sparsegrid(sg, f_on_grid, target_points)
    @test all(v[1] ≈ f(target_points[i]) for (i,v) in enumerate(interpolation_result))
end

@testset "Sparse Grid Integration Tests" begin
    f(x) = @. 3.0*x^2 + 2.0*x + 1.0
    pcl = precompute_lagrange_integrals(7)

    n,k = 2,0
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    f_on_grid = [[f(x[1]), f(x[2])^2] for x in get_grid_points(sg)]
    # Test integrate_on_sparsegrid
    integral_result = integrate_on_sparsegrid(sg, f_on_grid, pcl)
    @test integral_result ≈ [1.0,1.0]

    n,k = 2,4
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    f_on_grid = [[f(x[1]), f(x[2])^2] for x in get_grid_points(sg)]
    # Test integrate_on_sparsegrid
    integral_result = integrate_on_sparsegrid(sg, f_on_grid, pcl)

    @test integral_result ≈ [2.0,92/15]

    # Test integrate_L2_on_sparsegrid
    n,k = 2,5
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    f_on_grid = [[1.0] for x in get_grid_points(sg)]
    l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
    @test l2_integral_result[1] ≈ 1.0

    n,k = 1,1
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    f_on_grid = [[x[1]] for x in get_grid_points(sg)]
    pairwise_norms = precompute_pairwise_norms(f_on_grid, (x,y)->dot(x,y))
    @test pairwise_norms ≈ [0.0 0.0 0.0; 0.0 1.0 -1.0; 0.0 -1.0 1.0]

    l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
    @test l2_integral_result[1] ≈ sqrt(1/3)

    n,k = 1,3
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    f_on_grid = [[x[1]^2] for x in get_grid_points(sg)]
    l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
    @test l2_integral_result[1] ≈ sqrt(1/5)

    # Test precompute_pairwise_norms
    n,k = 2,4
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    f_on_grid = [[f(x[1])] for x in get_grid_points(sg)]
    pairwise_norms = precompute_pairwise_norms(f_on_grid, (x,y)->x.*y)
    l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
    @test l2_integral_result ≈ sqrt(92/15)
end
