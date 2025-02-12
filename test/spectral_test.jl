@testset "Spectral Representation Tests" begin
    @testset "Genz" begin
        n = 8
        C = 1.0
        W = 0.0
        T = "exponentialdecay"
        N = "gaussianpeak"
        f = genz(n::Int, C::Float64, W::Float64, T::String, N::String)

        # Approximate
        (sg, f_on_Z) = adaptive_sparsegrid(f, n)

        # Convert to spectral
        f_spectral = convert_to_spectral_approximation(sg, f_on_Z)

        sg_test = create_sparsegrid(create_smolyak_miset(n,4))
        test_points = get_grid_points(sg_test)

        f_spectral_on_test = f_spectral.(test_points)

        @test all(isapprox(f(x), f_spectral_on_test[i]; atol=1e-3) for (i,x) in enumerate(test_points))
    end
end