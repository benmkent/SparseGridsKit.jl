@testset "Sparse Grid Adaptivity Tests" begin
    @testset "One dimensional" begin
        test_points = range(-1,stop=1,length=100)

        f(x) = @. 3.0*x[1]^2 + 2.0*x[1] + 1.0
        nparams = 1

        (sg, f_on_Z) = adaptive_sparsegrid(f, nparams)
        # Expect three point approx cubic (1 iteration to suffice)
        @test all([f(x)] ≈ interpolate_on_sparsegrid(sg,f_on_Z,x) for x in test_points)
        @test get_n_grid_points(sg) == 3

        f2(x) = f(x).^2
        (sg, f2_on_Z) = adaptive_sparsegrid(f2, nparams)
        # Expect three point approx cubic (1 iteration to suffice)
        @test all([f2(x)] ≈ interpolate_on_sparsegrid(sg,f2_on_Z,x) for x in test_points)
        @test get_n_grid_points(sg) == 5

        f3(x) = f(x).^3
        (sg, f3_on_Z) = adaptive_sparsegrid(f3, nparams)
        # Expect three point approx cubic (1 iteration to suffice)
        @test all([f3(x)] ≈ interpolate_on_sparsegrid(sg,f3_on_Z,x) for x in test_points)
        @test get_n_grid_points(sg) == 9
    end
    @testset "Genz" begin
        n = 8
        C = 1.0
        W = 0.0
        T = "exponentialdecay"
        N = "gaussianpeak"
        f = genz(n::Int, C::Float64, W::Float64, T::String, N::String)

        # Approximate
        (sg, f_on_Z) = adaptive_sparsegrid(f, n)
        sg_test = create_sparsegrid(create_smolyak_miset(n,3))
        test_points = get_grid_points(sg_test)

        f_approx_on_test = interpolate_on_sparsegrid(sg,f_on_Z,test_points)
        
        @test all(isapprox(f(x), f_approx_on_test[i]; atol=1e-3) for (i,x) in enumerate(test_points))

        # Test termination by max points
        test_max_points = 10
        (sg, f_on_Z) = adaptive_sparsegrid(f, n, maxpts=test_max_points)
        @test isapprox(get_n_grid_points(sg), test_max_points, atol=10)
    end
end