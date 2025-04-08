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
    @test sort(grid_points) ≈ sort([   [0.0],
                            [1.0],
                            [-1.0],
                            [0.7071067811865475],
                            [-0.7071067811865475],
                            [0.9238795325112867],
                            [0.3826834323650898],
                            [-0.3826834323650898],
                            [-0.9238795325112867]])

    # Test map_from_to
    sg_from = create_sparsegrid(create_smolyak_miset(1, 8))
    sg_to = create_sparsegrid(create_smolyak_miset(1, 9))
    map_vector = map_from_to(sg_from, sg_to)
    @test get_grid_points(sg_from) ≈ get_grid_points(sg_to)[map_vector]
end