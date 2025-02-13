using ApproxFun
using SparseArrays
import Base.+

mutable struct SpectralSparseGridApproximation
    dims::Int
    expansiondimensions::Vector{Int}
    polytypes::Vector{Any}
    coefficients::SparseVector{Any}
    polydegrees::Vector{Vector}
end

function SpectralSparseGridApproximation(dims, expansiondims, polytypes, coefficients)
    coefficients = sparse(coefficients)
    polydegrees, coefficients = get_spectral_poly_representation(expansiondims, coefficients)
    return SpectralSparseGridApproximation(dims, expansiondims, polytypes, coefficients, polydegrees)
end

function (s::SpectralSparseGridApproximation)(x)
    @assert length(x) == s.dims

    # Evaluate
        evaluation = nothing
        for (i,idx) in enumerate(s.coefficients.nzind)
            c = s.coefficients[idx]
            # Product of polynomials
            if isnothing(evaluation)
                evaluation = c .* prod(Fun(s.polytypes[k],[zeros(s.polydegrees[i][k]);1])(x[k]) for k = 1:s.dims)
            else
                evaluation = evaluation + c .* prod(Fun(s.polytypes[k],[zeros(s.polydegrees[i][k]);1])(x[k]) for k = 1:s.dims)
            end
        end
    return evaluation
end

function +(grid1::SpectralSparseGridApproximation, grid2::SpectralSparseGridApproximation)
    # Ensure dimensions and grid sizes match
    if grid1.dims != grid2.dims
        throw(DimensionMismatch("Tensors have different numbers of dimensions"))
    end
    
    # Ensure the polynomial types match for each dimension
    for i in 1:grid1.dims
        if typeof(grid1.polytypes[i]) != typeof(grid2.polytypes[i])
            throw(DimensionMismatch("Polynomial types for dimension $i do not match"))
        end
    end

    # Find pairwise maximums of dimensions
    max_dims = max.(grid1.expansiondimensions, grid2.expansiondimensions)

    # Unpack and pad coefficients
    coeffs1 = grid1.coefficients
    coeffs2 = grid2.coefficients

    # Create sparse zero vectors of the product of the max dimensions
    coeffs1_padded = spzeros(Float64, prod(max_dims))
    coeffs2_padded = spzeros(Float64, prod(max_dims))


    # Map original coefficients to the new padded vector for grid1 (sparse array handling)
    if issparse(coeffs1)
        for i in 1:length(coeffs1.nzval)
            original_index = coeffs1.nzind[i]  # Get the sparse index
            indices = reverse_kron_index(original_index, grid1.expansiondimensions)  # Get multi-dimensional index
            new_index = kron_index(indices, max_dims)  # Convert multi-dimensional index to linear index
            coeffs1_padded[new_index] = coeffs1.nzval[i]  # Assign value at the correct index
        end
    else
        for i in 1:length(coeffs1)
            indices = reverse_kron_index(i, grid1.expansiondimensions)
            new_index = kron_index(indices, max_dims)
            coeffs1_padded[new_index] = coeffs1[i]
        end
    end

    # Map original coefficients to the new padded vector for grid2 (sparse array handling)
    if issparse(coeffs2)
        for i in 1:length(coeffs2.nzval)
            original_index = coeffs2.nzind[i]  # Get the sparse index
            indices = reverse_kron_index(original_index, grid2.expansiondimensions)  # Get multi-dimensional index
            new_index = kron_index(indices, max_dims)  # Convert multi-dimensional index to linear index
            coeffs2_padded[new_index] = coeffs2.nzval[i]  # Assign value at the correct index
        end
    else
        for i in 1:length(coeffs2)
            indices = reverse_kron_index(i, grid2.expansiondimensions)
            new_index = kron_index(indices, max_dims)
            coeffs2_padded[new_index] = coeffs2[i]
        end
    end

    # Return the sum of the two sparse vectors
    return SpectralSparseGridApproximation(grid1.dims, max_dims, grid1.polytypes, coeffs1_padded + coeffs2_padded)
end

function -(grid1::SpectralSparseGridApproximation, grid2::SpectralSparseGridApproximation)
    grid2.coefficients = -grid2.coefficients
    return grid1 + grid2 
end 

function get_spectral_poly_representation(expansiondimensions, coefficients)
        # Create an empty array to store the tensor index and corresponding coefficient
        poly = Vector{Vector{Int}}(undef,length(coefficients.nzval))

        # Iterate over non-zero coefficients
        for i in 1:length(coefficients.nzval)
            # Get the linear index of the non-zero coefficient
            linear_index = coefficients.nzind[i]
            
            # Convert the linear index to multi-dimensional indices using reverse_kron_index
            multi_dim_index = reverse_kron_index(linear_index, expansiondimensions)
                        
            # Store the tuple of the multi-dimensional index and coefficient
            poly[i] = multi_dim_index .- 1
        end
    
        return poly, coefficients
end

function poly_Fun(knots, ind; lb=-1.0, ub=1.0)
    n = length(knots)
    space = Chebyshev(lb..ub)

    f_knots = zeros(eltype(knots), length(knots))
    f_knots[ind] = 1

    V = Array{Float64}(undef,n,n);  # Create a Vandermonde matrix by evaluating the basis at the grid

    for k = 1:n
        V[:,k] = Fun(space,[zeros(k-1);1]).(knots)
    end

    f = Fun(space,V\f_knots);

    return f
end

function createlagrangepolys_approxfun(pts)
    p = Vector{Any}(undef, length(pts))
    for ii in eachindex(pts)
        p[ii] = poly_Fun(pts, ii)
    end
    return p
end

function convert_to_spectral_approximation(sparsegrid, fongrid)

    terms = sparsegrid.terms
    cterms = sparsegrid.cterms
    maprowtouniquept = sparsegrid.maptermstosparsegrid

    ndims = sparsegrid.dims
    
    if isa(sparsegrid.polyperlevel[1][1],Function)
        poly = Vector{Any}(undef,ndims) 
        for ii = 1:ndims
            poly[ii] = sparsegrid.polyperlevel
        end
    else
        poly = sparsegrid.polyperlevel
    end


    polytypes = Vector(undef,ndims)
    for ii = 1:ndims
        polytypes[ii] = poly[ii][1][1].space
    end

    SSG_total = SpectralSparseGridApproximation(ndims, zeros(Int,ndims), polytypes, [])
    for (i, row) = enumerate(eachrow(terms))
        if cterms[i] == 0
            continue
        end
        vectors  = [poly[j][row[1][j][1]][row[1][j][2]].coefficients for j = 1:ndims]
        f_Fun_i, dims = truncated_kron(vectors)
        SSG_i = SpectralSparseGridApproximation(ndims, dims, polytypes, fongrid[maprowtouniquept[i]] * cterms[i] * f_Fun_i)
        SSG_total = SSG_i + SSG_total
    end

    return SSG_total
end

function truncated_kron(vectors; tol=1e-10)
    result = sparse(vectors[1])
    for v in vectors[2:end]
        result = kron(result, sparse(v))

        # Apply truncation
        result = sparse(result .* (abs.(result) .> tol))
    end

    # Reshape into a tensor
    tensor_dims = length.(vectors)
    # return reshape(result, tensor_dims...)
    return result, tensor_dims
end

function kron_index(indices::Vector{Int}, dims::Vector{Int})
   # Convert to 0-based indices for easier calculation
    indices_0based = indices .- 1
    index = 0
    for ii = 1:length(indices)
        index = index + indices_0based[ii]*prod(dims[ii+1:end])
    end
    return index + 1
end

function reverse_kron_index(index::Int, dims::Vector{Int})
        index_0based = index - 1  # Convert to 0-based indexing
        indices = Vector{Int}(undef, length(dims))
        for ii = 1:length(dims)
            prod_rest = prod(dims[ii+1:end])  # Product of remaining dimensions
            indices[ii] = (index_0based ÷ prod_rest) + 1
            index_0based %= prod_rest
        end
        return indices   
end