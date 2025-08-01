using ApproxFun
using SparseArrays
using Polynomials
import Base.+

"""
    SpectralSparseGridApproximation

A spectral sparse grid approximation.

# Fields
- `dims::Int`: Dimension of the domain.
- `expansiondimensions::Vector{Int}`: Vector of maximum polynomial degree in each dimension.
- `polytypes::Vector{Any}`: Vector of polynomial spaces.
- `coefficients::SparseVector{Any}`: SparseVector of polynomial coefficients
- `polydegrees::SparseVector{Vector}`: SparseVector of polynomial degrees for each non zero coefficient.
"""
mutable struct SpectralSparseGridApproximation
    dims::Int
    expansiondimensions::Vector{Int}
    polytypes::Vector{Any}
    coefficients::SparseVector
    polydegrees
    domain
end

"""
    SpectralSparseGridApproximation(dims, expansionparams, polytypes, coefficients)

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
function SpectralSparseGridApproximation(dims, expansionparams, polytypes, coefficients)
    coefficients = sparse(coefficients)
    polydegrees, coefficients = get_spectral_poly_representation(expansionparams, coefficients)
    domain = Vector{Vector{Real}}(undef, dims)
    for i in 1:length(polytypes)
        domain[i] = getdomain(polytypes[i])
    end
    return SpectralSparseGridApproximation(dims, expansionparams, polytypes, coefficients, polydegrees,domain)
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
        evaluation = 0.0*s.coefficients[s.coefficients.nzind[1]]
        for (i,idx) in enumerate(s.coefficients.nzind)
            c = s.coefficients[idx]
            # Product of polynomials
            # if isnothing(evaluation)
                # evaluation = c .* prod(Fun(s.polytypes[k],[zeros(s.polydegrees[idx][k]);1])(x[k]) for k = 1:s.dims)
            # else
                # evaluation = evaluation + c .* prod(Fun(s.polytypes[k],[zeros(s.polydegrees[idx][k]);1])(x[k]) for k = 1:s.dims)
                evaluation = evaluation + c .* prod(Fun(s.polytypes[k],[zeros(s.polydegrees[idx][k]);1])(x[k]) for k = 1:s.dims)

            # end
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
    coeffs1_padded = spzeros(eltype(coeffs1), prod(max_dims))
    coeffs2_padded = spzeros(eltype(coeffs2), prod(max_dims))


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
    # Avoid having to create zeros.
    coeffs_combined = coeffs1_padded
    for idx in coeffs2_padded.nzind
        if idx ∈ coeffs1_padded.nzind
            coeffs_combined[idx] = coeffs1_padded[idx] + coeffs2_padded[idx]
        else
            coeffs_combined[idx] = coeffs2_padded[idx]
        end
    end
    return SpectralSparseGridApproximation(grid1.dims, max_dims, grid1.polytypes, coeffs_combined)
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

        poly = SparseVector(length(coefficients), coefficients.nzind, poly)    
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
    terms = sparsegrid.data.terms
    coeff_per_term = sparsegrid.data.coeff_per_term
    maprowtouniquept = sparsegrid.data.terms_to_grid_points
    domain = sparsegrid.domain

    nparams = sparsegrid.dims
    
    if isa(sparsegrid.data.point_sequences_per_dimension, Vector{AbstractVector{<:Real}})
        point_sequences_per_dimension = Vector{Vector{Vector{<:Real}}}(undef,nparams) 
        for ii = 1:nparams
            point_sequences_per_dimension[ii] = sparsegrid.data.point_sequences_per_dimension
        end
    else
        point_sequences_per_dimension = sparsegrid.data.point_sequences_per_dimension
    end

    polytypes = Vector{Space}(undef,nparams)
    for ii = 1:nparams
        # polytypes[ii] = poly[ii][1][1].space
        polytypes[ii] = Space(Interval(domain[ii]...))
    end

    SSG_total = nothing
    coeff_ij = Vector{Vector{Float64}}(undef, nparams)
    for (i, row) = enumerate(eachrow(terms))
        if coeff_per_term[i] == 0
            continue
        end
        for j = 1:nparams
            S = polytypes[j]
            p = points(S,length(point_sequences_per_dimension[j][end]))
            v = Vector{Float64}(undef, length(p))
            lagrange_evaluation!(v, point_sequences_per_dimension[j][row[1][j][1]],row[1][j][2],p);  # values at the default grid
            f = Fun(S,ApproxFun.transform(S,v));
            coeff_ij[j] = f.coefficients
        end
        # vectors  = [poly[j][row[1][j][1]][row[1][j][2]].coefficients for j = 1:nparams]
        f_Fun_i, dims = truncated_kron(coeff_ij; tol=1e-15)
        # vector_data = Vector{eltype(fongrid)}(undef, length(f_Fun_i.nzind))

        # for j in eachindex(f_Fun_i.nzind)
        #     # Get the multi-index for the current non-zero coefficient
        #     vector_data[j] = coeff_per_term[i] * fongrid[maprowtouniquept[i]] * f_Fun_i[j]
        # end
        sp_data = Vector{eltype(fongrid)}(undef,length(f_Fun_i.nzind))
        for (ii, sparse_ii) in enumerate(f_Fun_i.nzind)
            # Get the multi-index for the current non-zero coefficient
            sp_data[ii] = coeff_per_term[i] * fongrid[maprowtouniquept[i]] * f_Fun_i[sparse_ii]
        end
        # sp_data = coeff_per_term[i] * fongrid[maprowtouniquept[i]] * f_Fun_i;
        # SparseVector(f_Fun_i.n,f_Fun_i.nzind, coeff_per_term[i] * fongrid[maprowtouniquept[i]] * f_Fun_i)
        SSG_i = SpectralSparseGridApproximation(nparams, dims, polytypes, sparsevec(f_Fun_i.nzind,sp_data))
        if isnothing(SSG_total)
            SSG_total = SSG_i
        else
            # Add to the existing SSG_total
            SSG_total = SSG_total + SSG_i
        end
    end

    return SSG_total
end

"""
    convert_to_sg_approximation(ssg::SpectralSparseGridApproximation)
    Convert a `SpectralSparseGridApproximation` object to a SparseGridApproximation representation.

# Arguments
- `ssg::SpectralSparseGridApproximation`: Spectral sparse grid approximation.
# Returns
- `SparseGridApproximation`: Sparse grid approximation representation of `ssg`.
"""
function convert_to_sg_approximation(ssg::SpectralSparseGridApproximation, knots, rules)
    # First identify the necessary polynomial spaces
    polydegrees = ssg.polydegrees

    # Create inverse maps
    inverse_maps = [inverse_level_map(rules[i], ssg.expansiondimensions[i]) for i in 1:ssg.dims]

    # Next identify the required multi-index set
    multi_indices = Vector{Vector{Int}}(undef, length(ssg.polydegrees.nzind))
    for (i, idx) in enumerate(ssg.polydegrees.nzind)
        pd = polydegrees[idx]
        # Convert to multi-index
        multi_indices[i] = [inverse_maps[i][pd[i]+1] for i in 1:ssg.dims]
    end
    miset = MISet(multi_indices)
    adm, miss = check_admissibility(miset)
    miset_admissible = add_mi(miset,miss)

    # Evaluate the SpectralSparseGridApproximation at the knots of the sparse grid
    sg = create_sparsegrid(miset_admissible, rule=rules, knots=knots)

    # Create the SparseGridApproximation object
    eval_on_grid = ssg.(get_grid_points(sg))
    
    return SparseGridApproximation(sg, eval_on_grid)
end
"""
    inverse_level_map(levelfunction, maxlevel)
    Inverse level map for a given level2knots type function.

# Arguments
- `levelfunction`: Level function.
- `maxlevel`: Maximum level.
# Returns
- `Vector{Int}`: Inverse level map.
"""
function inverse_level_map(levelfunction, maxlevel)
    # Create a vector to store the inverse level map
    inverse_map = Vector{Int}(undef, levelfunction(maxlevel))
    # Iterate over the levels and fill the inverse map
    for ii = maxlevel:-1:1
        nknots = levelfunction(ii)
        for jj = 0:(nknots-1)
            inverse_map[jj+1] = ii
        end
    end
    return inverse_map
end
"""
    truncated_kron(vectors; tol=1e-10)

Compute the Kronecker product of a list of vectors with truncation.

# Arguments
- `vectors::Vector{Vector}`: Vector of vectors to compute the Kronecker product.
- `tol::Float64`: Tolerance value for coefficienet truncation to promote sparsity. Elements with absolute values less than `tol` are set to zero. Default is `1e-15` to remove coefficients at machine precision.

# Returns
- `result`: SparseVector representation of truncated Kronecker product of the input vectors.
- `tensor_dims::Vector{Int}`: Dimensions of the resulting tensor.
"""
function truncated_kron(vectors; tol=1e-15)
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

"""
    getdomain(polytype)

    Takes an ApproxFun polynomial space a returns a vector [a,b] of the domain endpoints.

# Arguments
- `polytype` : ApproxFun Space

# Returns
- `domain` : Vector domain [a,b]
"""
function getdomain(polytype)
    d = domain(polytype)
    a,b = leftendpoint(d), rightendpoint(d)        
    return [a,b]
end