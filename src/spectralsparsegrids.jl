using ApproxFun
using SparseArrays
import Base.+

"""
    SpectralSparseGridApproximation

A spectral sparse grid approximation.

# Fields
- `dims::Int`: Dimension of the domain.
- `expansiondimensions::Vector{Int}`: Vector of maximum polynomial degree in each dimension.
- `polytypes::Vector{Any}`: Vector of polynomial spaces.
- `coefficients::SparseVector{Any}`: Sparse ector of polynomial coefficients
- `polydegrees::SparseVector{Vector}`: SparseVector of polynomial degrees for each non zero coefficient.
"""
mutable struct SpectralSparseGridApproximation
    dims::Int
    expansiondimensions::Vector{Int}
    polytypes::Vector{Any}
    coefficients::SparseVector{Any}
    polydegrees
    domain
end

"""
    SpectralSparseGridApproximation(dims, expansiondims, polytypes, coefficients)

Constructs a spectral sparse grid approximation.

# Arguments
- `dims::Int`: Dimension of the domain.
- `expansiondimensions::Vector{Int}`: Vector of maximum polynomial degree in each dimension.
- `polytypes::Vector{Any}`: Vector of polynomial spaces.
- `coefficients::Vector{Any}`: Kronecker product vector of polynomial coefficients
- `domain::Vector{Vector{Real}}`: Vector of domains for inputs

# Returns
- `SpectralSparseGridApproximation`: An instance of `SpectralSparseGridApproximation`.
"""
function SpectralSparseGridApproximation(dims, expansiondims, polytypes, coefficients, domain)
    coefficients = sparse(coefficients)
    polydegrees, coefficients = get_spectral_poly_representation(expansiondims, coefficients)
    return SpectralSparseGridApproximation(dims, expansiondims, polytypes, coefficients, polydegrees, domain)
end

"""
    (s::SpectralSparseGridApproximation)(x)

Evaluate the spectral sparse grid approximation at the given point `x`.

# Arguments
- `s::SpectralSparseGridApproximation`: Spectral sparse grid approximation object.
- `x::Vector{T}`: Vector x in domain to be evaluated

# Returns
- `evaluation`: Evaluation of s(x)
"""
function (s::SpectralSparseGridApproximation)(x)
    @assert length(x) == s.dims
    # Evaluate
        evaluation = nothing
        for (i,idx) in enumerate(s.coefficients.nzind)
            c = s.coefficients[idx]
            # Product of polynomials
            if isnothing(evaluation)
                evaluation = c .* prod(Fun(s.polytypes[k],[zeros(s.polydegrees[idx][k]);1])(x[k]) for k = 1:s.dims)
            else
                evaluation = evaluation + c .* prod(Fun(s.polytypes[k],[zeros(s.polydegrees[idx][k]);1])(x[k]) for k = 1:s.dims)
            end
        end
    return evaluation
end

"""
    +(grid1::SpectralSparseGridApproximation, grid2::SpectralSparseGridApproximation)

Add two `SpectralSparseGridApproximation` objects.

# Arguments
- `grid1::SpectralSparseGridApproximation`: Spectral sparse grid approximation.
- `grid2::SpectralSparseGridApproximation`: Second spectral sparse grid approximation.

# Returns
- `SpectralSparseGridApproximation` representing sum of `grid1` and `grid2`.
"""
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
    return SpectralSparseGridApproximation(grid1.dims, max_dims, grid1.polytypes, coeffs1_padded + coeffs2_padded, grid1.domain)
end

"""
    -(grid1::SpectralSparseGridApproximation, grid2::SpectralSparseGridApproximation)

Subtracts `SpectralSparseGridApproximation` objects

# Arguments
- `grid1::SpectralSparseGridApproximation`: Spectral sparse grid approximation.
- `grid2::SpectralSparseGridApproximation`: Second spectral sparse grid approximation, to be subtracted

# Returns
- `SpectralSparseGridApproximation` object representing the subtraction.
"""
function -(grid1::SpectralSparseGridApproximation, grid2::SpectralSparseGridApproximation)
    grid2copy = deepcopy(grid2)
    grid2copy.coefficients = -grid2copy.coefficients
    return grid1 + grid2copy 
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

        poly = sparse(coefficients.nzind, ones(length(coefficients.nzind)), poly, length(coefficients), 1)
    
        return poly, coefficients
end

"""
    poly_Fun(knots, ind; lb=-1.0, ub=1.0)

Lagrange interpolation polynomial `ind` for the set of `knots`

# Arguments
- `knots::Vector{T}`: A vector of knots
- `ind::Int`: Knot index for polynomial equal to 1.
- `lb::Float64`: Domain lower bound. Default is -1.0.
- `ub::Float64`: Domain upper bound. Default is 1.0.

# Returns
- `f::Fun`: Lagrange interpolation polynomial as ApproxFun Fun object.
"""
function poly_Fun(knots, ind; space=Chebyshev(-1.0..1.0))
    n = length(knots)

    f_knots = zeros(eltype(knots), length(knots))
    f_knots[ind] = 1

    V = Array{Float64}(undef,n,n);  # Create a Vandermonde matrix by evaluating the basis at the grid

    for k = 1:n
        V[:,k] = Fun(space,[zeros(k-1);1]).(knots)
    end

    f = Fun(space,V\f_knots);

    return f
end

"""
    createlagrangepolys_approxfun(pts)

Create a vector of Lagrange polynomial approximation functions for a given set of knots.

# Arguments
- `pts`: Vector of knots.

# Returns
- `Vector{Any}`: Vector containing the Lagrange polynomial approximation `Fun`s for each point in `pts`.

"""
function createlagrangepolys_approxfun(pts, space=Chebyshev(-1.0..1.0))
    p = Vector{Any}(undef, length(pts))
    for ii in eachindex(pts)
        p[ii] = poly_Fun(pts, ii, space)
    end
    return p
end

"""
    convert_to_spectral_approximation(sga::SparseGridApproximation)

Convert a `SparseGridApproximation` object to its spectral approximation.

# Arguments
- `sga::SparseGridApproximation`: Sparse grid approximation.

# Returns
- Spectral Sparse Grid Approximation representing `sga`.
"""
function convert_to_spectral_approximation(sga::SparseGridApproximation)
    return convert_to_spectral_approximation(sga.sg, sga.fongrid)
end

"""
    convert_to_spectral_approximation(sparsegrid::SparseGrid, fongrid)

Converts a `sparse grid` and corresponding evaluations `fongrid` to a SpectralSparseGridApproximation

# Arguments
- `sparsegrid::SparseGrid`: Sparse grid
- `fongrid`: Function evaluations of `sparse grid`.

# Returns
- `SpectralSparseGridApproximation` represntation of `sparsegrid` with evaluations `fongrid`.
"""
function convert_to_spectral_approximation(sparsegrid::SparseGrid, fongrid)
    terms = sparsegrid.terms
    cterms = sparsegrid.cterms
    maprowtouniquept = sparsegrid.maptermstosparsegrid
    domain = sparsegrid.domain

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

    SSG_total = SpectralSparseGridApproximation(ndims, zeros(Int,ndims), polytypes, [], domain)
    for (i, row) = enumerate(eachrow(terms))
        if cterms[i] == 0
            continue
        end
        vectors  = [poly[j][row[1][j][1]][row[1][j][2]].coefficients for j = 1:ndims]
        f_Fun_i, dims = truncated_kron(vectors)
        SSG_i = SpectralSparseGridApproximation(ndims, dims, polytypes, fongrid[maprowtouniquept[i]] * cterms[i] * f_Fun_i, domain)
        SSG_total = SSG_i + SSG_total
    end

    return SSG_total
end

"""
    truncated_kron(vectors; tol=1e-10)

Compute the Kronecker product of a list of vectors with truncation.

# Arguments
- `vectors::Vector{Vector}`: Vector of vectors to compute the Kronecker product.
- `tol::Float64`: Tolerance value for coefficienet truncation to promote sparsity. Elements with absolute values less than `tol` are set to zero. Default is `1e-10`.

# Returns
- `result`: SparseVector representation of truncated Kronecker product of the input vectors.
- `tensor_dims::Vector{Int}`: Dimensions of the resulting tensor.
"""
function truncated_kron(vectors; tol=1e-10)
    result = sparse(vectors[1])
    for v in vectors[2:end]
        result = kron(result, sparse(v))

        # Apply truncation
        result = sparse(result .* (abs.(result) .> tol))
    end

    tensor_dims = length.(vectors)
    return result, tensor_dims
end

"""
    kron_index(indices::Vector{Int}, dims::Vector{Int})

Calculate the linear index in a Kronecker product space given the multi-dimensional indices and the dimensions of each space.

# Arguments
- `indices::Vector{Int}`: Index vector.
- `dims::Vector{Int}`: Vector of storage dimension for each tensor dimension.

# Returns
- Linear index corresponding to the index given in `indices`.

"""
function kron_index(indices::Vector{Int}, dims::Vector{Int})
   # Convert to 0-based indices for easier calculation
    indices_0based = indices .- 1
    index = 0
    for ii = 1:length(indices)
        index = index + indices_0based[ii]*prod(dims[ii+1:end])
    end
    return index + 1
end

"""
    reverse_kron_index(index::Int, dims::Vector{Int})

Given linear index `index` and a vector of dimensions `dims`, computes the multi-dimensional tensor index.

# Arguments
- `index::Int`: Linear index.
- `dims::Vector{Int}`: Vector containing the storage dimension of the tensor dimensions.

# Returns
- `Vector{Int}`: Tensor index

"""
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