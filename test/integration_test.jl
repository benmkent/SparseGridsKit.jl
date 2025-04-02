@testset "Sparse Grid Integration Tests" begin
    f(x) = @. 3.0*x^2 + 2.0*x + 1.0

    n,k = 2,0

    domain=fill([-1,1],n)
    pcl = precompute_lagrange_integrals(7, domain, CCPoints(), Doubling())

    mi_set = create_smolyak_miset(n,k)
    domain = fill([-1, 1], n)
    sg = create_sparsegrid(mi_set,domain)
    f_on_grid = [[f(x[1]), f(x[2])^2] for x in get_grid_points(sg)]
    # Test integrate_on_sparsegrid
    integral_result = integrate_on_sparsegrid(sg, f_on_grid, pcl)
    @test integral_result ≈ [1.0,1.0]

    n,k = 2,4
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set,domain)
    f_on_grid = [[f(x[1]), f(x[2])^2] for x in get_grid_points(sg)]
    # Test integrate_on_sparsegrid
    integral_result = integrate_on_sparsegrid(sg, f_on_grid, pcl)

    @test integral_result ≈ [2.0,92/15]

    # Test integrate_L2_on_sparsegrid
    n,k = 2,5
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set,domain)
    f_on_grid = [[1.0] for x in get_grid_points(sg)]
    l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
    @test l2_integral_result[1] ≈ 1.0

    n,k = 1,1
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set,domain)
    f_on_grid = [[x[1]] for x in get_grid_points(sg)]
    pairwise_norms = precompute_pairwise_norms(f_on_grid, product=(x,y)->dot(x,y))
    @test pairwise_norms ≈ [1.0 0.0 -1.0; 0.0 0.0 0.0; -1.0 0.0 1.0]

    f_on_grid = [x[1] for x in get_grid_points(sg)]
    pairwise_norms = precompute_pairwise_norms(f_on_grid)
    @test pairwise_norms ≈ [1.0 0.0 -1.0; 0.0 0.0 0.0; -1.0 0.0 1.0]

    l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
    @test l2_integral_result ≈ sqrt(1/3)

    n,k = 1,3
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set,domain)
    f_on_grid = [x[1]^2 for x in get_grid_points(sg)]
    l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
    @test l2_integral_result ≈ sqrt(1/5)

    # Test precompute_pairwise_norms and vector
    n,k = 2,4
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set,domain)
    f_on_grid = [[f(x[1])] for x in get_grid_points(sg)]
    pairwise_norms = precompute_pairwise_norms(f_on_grid)
    l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
    @test l2_integral_result ≈ [sqrt(92/15)]
end