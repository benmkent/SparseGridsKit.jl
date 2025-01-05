@testset "Sparse Grid Interpolation Test" begin
    # Test interpolate_on_sparsegrid
    f(x) = @. 3.0*x[1]^2 + 2.0*x[1] + 1.0

    n,k = 2,0
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    f_on_grid = [[f(x), cos(x[2])] for x in get_grid_points(sg)]
    target_points = [[x, x] for x in -1.0:0.1:1.0]
    interpolation_result = interpolate_on_sparsegrid(sg, f_on_grid, target_points)
    @test all(v ≈ [1.0, 1.0] for v in interpolation_result)

    n,k = 2,1
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    f_on_grid = [[f(x), cos(x[2])] for x in get_grid_points(sg)]
    target_points = [[x, x] for x in -1.0:0.1:1.0]
    interpolation_result = interpolate_on_sparsegrid(sg, f_on_grid, target_points)
    @test all(v[1] ≈ f(target_points[i]) for (i,v) in enumerate(interpolation_result))

    n,k = 4,3
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    # Define a complicated function
    ndims = 400
    f(x) = real.([(2.0 .+ (cos(2*i/ndims*prod(1.0.+x)))) for i=1:ndims])
    f_on_grid = [f(x) for x in get_grid_points(sg)]
    mi_enriched = add_mi(mi_set, get_reduced_margin(mi_set))
    sg_enriched = create_sparsegrid(mi_enriched)
    f_on_grid_2 = interpolate_on_sparsegrid(sg,f_on_grid,get_grid_points(sg_enriched))
    # Both should represent the same polynomial
    nmc = Integer(1e2);
    v = [2*(rand(n).-1) for ii=1:nmc]
    f_on_v = interpolate_on_sparsegrid(sg,f_on_grid,v)
    f2_on_v = interpolate_on_sparsegrid(sg_enriched,f_on_grid_2,v)
    @test all(isapprox(f_on_v[ii],f2_on_v[ii],rtol=1e-6) for ii in eachindex(v))
end