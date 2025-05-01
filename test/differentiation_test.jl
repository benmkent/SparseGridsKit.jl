@testset "Differentiation" begin
    # Test derivative function   
    # Create a sparse grid approximation of a polynomial
    n = 1
    p(x) = 5*x[1].^5 + x[1].^3 + 2*x[1].^2 + x[1] + 1
    p_prime(x) = 25*x[1].^4 + 3*x[1].^2 + 4*x[1] + 1

    sg = create_sparsegrid(create_smolyak_miset(n, 4))
    f_on_Z = p.(get_grid_points(sg))
    f_sga = SparseGridApproximation(sg, f_on_Z)
    f_sg_diff = derivative(f_sga)

    # Create test points
    test_points = [[x] for x in range(-1, stop=1, length=100)]
    # Evaluate the original function and its derivative at the test points
    f_test = p_prime.(test_points)
    f_test_diff = f_sg_diff.(test_points)

    @test all(isapprox(f_test[i], f_test_diff[i][1]) for i in 1:length(test_points))

    # Test alternative dispathc
    sg, f_sg_diff = derivative(sg, f_on_Z)
    f_test_diff = f_sg_diff.(test_points)
    @test all(isapprox(f_test[i], f_test_diff[i][1]) for i in 1:length(test_points))    
    
    # Test Spectral dispatch
    ssga = convert_to_spectral_approximation(f_sga)
    ssga_diff = derivative(ssga)
    f_test_diff = ssga_diff.(test_points)
    @test all(isapprox(f_test[i], f_test_diff[i][1]) for i in 1:length(test_points))

    # Test multi-variate
    n = 2
    p(x) = [3*x[1]^3*x[2]^2]
    p_prime(x) = [9*x[1]^2*x[2]^2, 6*x[1]^3*x[2]]

    sg = create_sparsegrid(create_smolyak_miset(n, 4))
    f_on_Z = p.(get_grid_points(sg))
    f_sga = SparseGridApproximation(sg, f_on_Z)
    f_sg_diff = derivative(f_sga)

    # Create test points
    test_points = [[x,y] for x in range(-1, stop=1, length=10) for y in range(-1, stop=1, length=10)]
    # Evaluate the original function and its derivative at the test points
    f_test = p_prime.(test_points)
    f_test_diff = f_sg_diff.(test_points)

    @test all(isapprox(f_test[i], f_test_diff[i]) for i in 1:length(test_points))
end

