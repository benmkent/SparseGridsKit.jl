using SparseGridsKit: truncated_kron, kron_index, reverse_kron_index
using SparseArrays: sparse
using ApproxFun

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

    @testset "Addition" begin
        space = Chebyshev(-1..1)
        tensor_1, dims_1 = truncated_kron([[1.0,2.0,3.0],[1.0]])
        tensor_2, dims_2 = truncated_kron([[1.0],[1.0,2.0,3.0]])
        domain = fill([-1,1],2)
        ssg1 = SpectralSparseGridApproximation(2,dims_1,[space,space],tensor_1,domain)
        ssg2 = SpectralSparseGridApproximation(2,dims_2,[space,space],tensor_2,domain)

        ssg = ssg1 + ssg2

        @test ssg.expansiondimensions == [3,3]
        @test ssg.polytypes == [space,space]
        @test all([ ssg.coefficients[1]  ==  2.0
                    ssg.coefficients[2]  ==  2.0
                    ssg.coefficients[3]  ==  3.0
                    ssg.coefficients[4]  ==  2.0
                    ssg.coefficients[7]  ==  3.0])
        @test all([ ssg.polydegrees[1] .== [0, 0]
                    ssg.polydegrees[2] .== [0, 1]
                    ssg.polydegrees[3] .== [0, 2]
                    ssg.polydegrees[4] .== [1, 0]
                    ssg.polydegrees[7] .== [2, 0]])

        
        ssg = ssg1 - ssg2
        # No [0,0] term after subtraction
        @test ssg.expansiondimensions == [3,3]
        @test ssg.polytypes == [space,space]
        @test all([ ssg.coefficients[2]  ==  -2.0
                    ssg.coefficients[3]  ==  -3.0
                    ssg.coefficients[4]  ==  2.0
                    ssg.coefficients[7]  ==  3.0])
        @test all([ ssg.polydegrees[2] .== [0, 1]
                    ssg.polydegrees[3] .== [0, 2]
                    ssg.polydegrees[4] .== [1, 0]
                    ssg.polydegrees[7] .== [2, 0]])
    end
    
    @testset "Representation" begin
        # Create simple SG
        n = 1
        C = 1.0
        W = 0.0
        T = "exponentialdecay"
        N = "gaussianpeak"
        f = genz(n::Int, C::Float64, W::Float64, T::String, N::String)
        domain = fill([-1,1],n)
        # Approximate
        (sg, f_on_Z) = adaptive_sparsegrid(f, domain, n)
        
        f_sg = SparseGridApproximation(sg,f_on_Z)

        # Convert to spectral
        f_spectral = convert_to_spectral_approximation(sg, f_on_Z)

        x_test = -1:0.1:1
        y_test = f_sg.(x_test)
        y_spectral = f_spectral.(x_test)
        @test all(isapprox(y_test[i], y_spectral[i]; atol=1e-8) for i in 1:length(x_test))

        n = 2
        domain = fill([-1,1],n)
        f = genz(n::Int, C::Float64, W::Float64, T::String, N::String)
        # Approximate
        (sg, f_on_Z) = adaptive_sparsegrid(f, domain, n)
        f_sg = SparseGridApproximation(sg,f_on_Z)

        # Convert to spectral
        f_spectral = convert_to_spectral_approximation(sg, f_on_Z)

        x_test = [rand(n) for i in 1:100]
        y_test = f_sg.(x_test)
        y_spectral = f_spectral.(x_test)
        @test all(isapprox(y_test[i], y_spectral[i]; atol=1e-8) for i in 1:length(x_test))

        n = 3
        f = genz(n::Int, C::Float64, W::Float64, T::String, N::String)
        domain = fill([-1,1],n)
        # Approximate
        (sg, f_on_Z) = adaptive_sparsegrid(f, domain, n)
        f_sg = SparseGridApproximation(sg,f_on_Z)

        # Convert to spectral
        f_spectral = convert_to_spectral_approximation(sg, f_on_Z)

        x_test = [rand(n) for i in 1:100]
        y_test = f_sg.(x_test)
        y_spectral = f_spectral.(x_test)
        @test all(isapprox(y_test[i], y_spectral[i]; atol=1e-8) for i in 1:length(x_test))

        n = 4
        f = genz(n::Int, C::Float64, W::Float64, T::String, N::String)
        domain = fill([-1,1],n)

        # Approximate
        (sg, f_on_Z) = adaptive_sparsegrid(f, domain, n)
        f_sg = SparseGridApproximation(sg,f_on_Z)

        # Convert to spectral
        f_spectral = convert_to_spectral_approximation(sg, f_on_Z)

        x_test = [rand(n) for i in 1:100]
        y_test = f_sg.(x_test)
        y_spectral = f_spectral.(x_test)
        @test all(isapprox(y_test[i], y_spectral[i]; atol=1e-8) for i in 1:length(x_test))

        n = 8
        f = genz(n::Int, C::Float64, W::Float64, T::String, N::String)
        domain = fill([-1,1],n)

        # Approximate
        (sg, f_on_Z) = adaptive_sparsegrid(f, domain, n)
        f_sg = SparseGridApproximation(sg,f_on_Z)

        # Convert to spectral
        f_spectral = convert_to_spectral_approximation(sg, f_on_Z)

        x_test = [rand(n) for i in 1:100]
        y_test = f_sg.(x_test)
        y_spectral = f_spectral.(x_test)
        @test all(isapprox(y_test, y_spectral; atol=1e-8))
    end
end