using SparseGridsKit: truncated_kron, kron_index, reverse_kron_index

@testset "Spectral Representation Tests" begin
    @testset "Indexing" begin
        n = 4
        vectors = [1:n,1:n,1:n]
        v_tensor, dims = truncated_kron(vectors)
        @test dims == [4,4,4]
        @test size(v_tensor) == (64,)

        @test all(v_tensor[kron_index([i,j,k],dims)] == i*j*k for i in 1:n for j in 1:n for k in 1:n)

        @test all(v_tensor[i] == prod(reverse_kron_index(i,dims)) for i in eachindex(v_tensor))

        vectors = [1:2,1:5,1:3]
        v_tensor, dims = truncated_kron(vectors)
        @test dims == [2,5,3]
        @test all(v_tensor[kron_index([i,j,k],dims)] == i*j*k for i in 1:2 for j in 1:5 for k in 1:3)

        @test all(v_tensor[i] == prod(reverse_kron_index(i,dims)) for i in eachindex(v_tensor))
    end
    
    @testset "Genz" begin
        # Create simple SG
        n = 1
        C = 1.0
        W = 0.0
        T = "exponentialdecay"
        N = "gaussianpeak"
        f = genz(n::Int, C::Float64, W::Float64, T::String, N::String)
        # Approximate
        (sg, f_on_Z) = adaptive_sparsegrid(f, n)
        
        f_sg = SparseGridApproximation(sg,f_on_Z)

        # Convert to spectral
        f_spectral = convert_to_spectral_approximation(sg, f_on_Z)

        x_test = -1:0.1:1
        y_test = f_sg.(x_test)
        y_spectral = f_spectral.(x_test)
        @test all(isapprox(y_test[i], y_spectral[i]; atol=1e-8) for i in 1:length(x_test))

        n = 2
        f = genz(n::Int, C::Float64, W::Float64, T::String, N::String)
        # Approximate
        (sg, f_on_Z) = adaptive_sparsegrid(f, n)
        f_sg = SparseGridApproximation(sg,f_on_Z)

        # Convert to spectral
        f_spectral = convert_to_spectral_approximation(sg, f_on_Z)

        x_test = [rand(n) for i in 1:100]
        y_test = f_sg.(x_test)
        y_spectral = f_spectral.(x_test)
        @test all(isapprox(y_test[i], y_spectral[i]; atol=1e-8) for i in 1:length(x_test))

        n = 3
        f = genz(n::Int, C::Float64, W::Float64, T::String, N::String)
        # Approximate
        (sg, f_on_Z) = adaptive_sparsegrid(f, n)
        f_sg = SparseGridApproximation(sg,f_on_Z)

        # Convert to spectral
        f_spectral = convert_to_spectral_approximation(sg, f_on_Z)

        x_test = [rand(n) for i in 1:100]
        y_test = f_sg.(x_test)
        y_spectral = f_spectral.(x_test)
        @test all(isapprox(y_test[i], y_spectral[i]; atol=1e-8) for i in 1:length(x_test))

        n = 4
        f = genz(n::Int, C::Float64, W::Float64, T::String, N::String)
        # Approximate
        (sg, f_on_Z) = adaptive_sparsegrid(f, n)
        f_sg = SparseGridApproximation(sg,f_on_Z)

        # Convert to spectral
        f_spectral = convert_to_spectral_approximation(sg, f_on_Z)

        x_test = [rand(n) for i in 1:100]
        y_test = f_sg.(x_test)
        y_spectral = f_spectral.(x_test)
        @test all(isapprox(y_test[i], y_spectral[i]; atol=1e-8) for i in 1:length(x_test))

        n = 8
        f = genz(n::Int, C::Float64, W::Float64, T::String, N::String)
        # Approximate
        (sg, f_on_Z) = adaptive_sparsegrid(f, n)
        f_sg = SparseGridApproximation(sg,f_on_Z)

        # Convert to spectral
        f_spectral = convert_to_spectral_approximation(sg, f_on_Z)

        x_test = [rand(n) for i in 1:100]
        y_test = f_sg.(x_test)
        y_spectral = f_spectral.(x_test)
        @test all(isapprox(y_test, y_spectral; atol=1e-8))
    end
end