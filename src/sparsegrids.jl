using LinearAlgebra
using StaticArrays
using Polynomials

import Base.:+
import Base.:-

"""
    SparseTerm{N}

Structure representing a sparse grid `term` with `N` dimensions.

"""
struct SparseTerm{N}
    data::NTuple{N, SVector{2,Int}} # Tuple containing `N` static vectors of two integers each.
end

"""
    sparse_grid_data{N}

Underlying data for a sparse grid of dimension `N`.
"""
mutable struct sparse_grid_data{N}
    terms::Vector{SparseTerm{N}}
    points_to_unique_indices::Vector{NTuple{N,Int}} # Mapping from terms to unique points
    terms_to_grid_points::Vector{Int} # Mapping from terms to grid points
    coeff_per_term::Vector{Int} # Combination coefficient per term
    terms_to_unique_points_per_dimension::Vector{Dict{Vector{Int},Int}} # Mapping from terms to unique points per dimension
    point_sequences_per_dimension::Vector{Vector{Vector{Float64}}}
    unique_points_per_dimension::Vector{Vector{Float64}} # Unique points per dimension
    weight_sequences_per_dimension::Vector{Vector{Vector{Float64}}}
end

"""
    SparseGrid

Sparse grid structure.
"""
mutable struct SparseGrid
    dims::Int # Number of dimensions
    domain::Vector{Vector{Float64}} # Grid domain
    multi_index_set::Matrix{Int} # Matrix of multi-index columns [β_1 β_2 ...]
    combination_coeff::Vector{Int} # Combination technique coefficients
    grid_points::Matrix{Real} # Sparse grid
    quadrature_weights::Vector{Float64} # Quadrature weights
    knots::Vector{Points}
    rules::Vector{Level}
    data::sparse_grid_data # Structure for underlying construction data
end

"""
    SparseGrid(dims, domain, multi_index_set, combination_coeff, grid_points, knots, rules, data)

Constructs a `SparseGrid` object with quadrature weights.
"""
function SparseGrid(dims, domain, multi_index_set, combination_coeff, grid_points, knots, rules, data)
    sg = SparseGrid(dims, domain, multi_index_set, combination_coeff, grid_points, [], knots, rules, data)
    compute_quadrature_weights!(sg)
    return sg
end

"""
    SparseGridApproximation

Function approximation consisting of sparse grid and function values on the grid.
"""
mutable struct SparseGridApproximation
    sg::SparseGrid
    fongrid::Vector
    domain::Vector
end

"""
    SparseGridApproximation

Constructs a SparseGridApproximation.
"""
function SparseGridApproximation(sg, fongrid)
    SparseGridApproximation(sg,fongrid,sg.domain)
end

"""
    (sga::SparseGridApproximation)(x)

Evaluate SparseGridApproximation at point `x`.
"""
function (sga::SparseGridApproximation)(x)
    @assert length(x) == sga.sg.dims
    return interpolate_on_sparsegrid(sga.sg, sga.fongrid, [x])[1]
end

"""
    Base.:+(sg1::SparseGrid, sg2::SparseGrid)

Add two `SparseGrid` objects by adding their common combination coefficients and merging their multi-index sets.
"""
function Base.:+(sg1::SparseGrid, sg2::SparseGrid)
    dims = sg1.dims
    @assert dims == sg2.dims

    sg1maxmi = length(sg1.data.point_sequences_per_dimension)
    sg2maxmi = length(sg2.data.point_sequences_per_dimension)

    if sg1maxmi >= sg2maxmi
        sg = copy(sg1)
        sga = copy(sg1)
        sgb = copy(sg2)
    else
        sg = copy(sg2)
        sga = copy(sg2)
        sgb = copy(sg1)
    end

    for (imi, mi) in enumerate(eachcol(sgb.multi_index_set))
        index = findfirst(mi == mi2 for mi2 in eachcol(sga.multi_index_set))
        if index !== nothing
            sg.combination_coeff[index] = sg.combination_coeff[index] + sgb.combination_coeff[imi]
        else
            push!(sg.combination_coeff, sgb.combination_coeff[imi])
            sg.multi_index_set = hcat(sg.multi_index_set, mi)
        end
    end

    (multi_index_set,combination_coeff) = sortmi(sg.multi_index_set, sg.combination_coeff)

    (terms, coeff_per_term, combination_coeff, multi_index_set) = sparsegridterms(multi_index_set, sg.data.point_sequences_per_dimension, combination_coeff=sg.combination_coeff)

    (points_to_unique_indices, grid_points, terms_to_grid_points) = sparsegridpoints(terms, sg.data.terms_to_unique_points_per_dimension, sg.data.unique_points_per_dimension, size(sg.multi_index_set, 1))

    data = sparse_grid_data{dims}(terms, points_to_unique_indices, terms_to_grid_points, coeff_per_term, terms_to_unique_points_per_dimension, sg.data.point_sequences_per_dimension, sg.data.unique_points_per_dimension, weight_sequences_per_dimension)

    sparsegrid = SparseGrid(dims, sg.domain, multi_index_set, combination_coeff, grid_points, sg.knots, sg.rules, data)
    return sparsegrid
end


"""
    -(sg1::SparseGrid, sg2::SparseGrid)

Subtracts two `SparseGrid` objects by negating `sg2` and adding to `sg1`.
"""
function Base.:-(sg1::SparseGrid, sg2::SparseGrid)
    sg2minus = copy(sg2)
    sg2minus.data.coeff_per_term = -sg2minus.data.coeff_per_term
    sg2minus.combination_coeff = -sg2minus.combination_coeff
    sg = copy(sg1) + sg2minus
    return sg
end

function Base.:copy(sg::SparseGrid)
    sgcopy = deepcopy(sg)
    return sgcopy
end

"""
    sparsegridprecompute(maxmi, knots=CCPoints(), rule=Doubling())

Precomputes weighted L^2 product integrals for each pair of Lagrange interpolation polynomials.

# Arguments
- `maxmi`: maxmium  level to consider
- `knots`: knot sequence for generating Lagrange polynomials
- `rule`: level-to-knots rule for generating Lagrange polynomials

# return
- `productintegrals`: Vector fo each dimension holding pcl[ii,jj,kk,ll] equal to weighted integral of
                      polynomial kk for level ii against polynomial ll for level jj.
"""
function sparsegridprecompute(maxmi, knots=CCPoints(), rule=Doubling())
    if isa(knots, Points)
        nparams=1
        knots = fill(knots, nparams)
    end
    if isa(rule, Level)
        rule = fill(rule, nparams)
    end
    productintegrals = Vector{Any}(undef, length(knots))
    for ii = eachindex(knots)
        (unique_points_per_dimension, point_sequences_per_dimension, terms_to_unique_points_per_dimension, weight_sequences_per_dimension) = sparsegridonedimdata(maxmi, knots[ii], rule[ii], knots[ii].domain)
        productintegrals[ii] = sparsegridproductintegrals(point_sequences_per_dimension, maxmi, knots[ii], rule[ii])
    end
    @debug "Computed weighted L^2 product integrals in each dimension"
    return (; productintegrals)
end

"""
    sparsegridproductintegrals(point_sequences_per_dimension, maxmi, knots, rule)

Computes integrals of pairwise polynomials.
"""
function sparsegridproductintegrals(point_sequences_per_dimension, maxmi, knots, rule)
    @debug "Computing product integrals"

    # Check uniqueness of points
    for pts in point_sequences_per_dimension
        check_unique_pts(pts)
    end

    # Allocate arrays
    productintegrals = zeros(Float64, maxmi, maxmi, length(point_sequences_per_dimension[end]), length(point_sequences_per_dimension[end]))

    x,w = knots(2*rule(maxmi))
    wp1x = similar(x)

    pevalx = zeros(Float64, maxmi, length(point_sequences_per_dimension[end]), length(x))
    for ii = 1:maxmi
        for jj = 1:length(point_sequences_per_dimension[ii])
            write_location = @view pevalx[ii, jj, :]
            pts_input = point_sequences_per_dimension[ii]
            lagrange_evaluation!(write_location, pts_input, jj, x)
        end
    end

    for ii = 1:maxmi
        # Compute integrals
        for ll = 1:maxmi
            @inbounds for p1 in eachindex(point_sequences_per_dimension[ii])
                wp1x .= w .* pevalx[ii, p1, :]
                @inbounds for p2 in eachindex(point_sequences_per_dimension[ll])
                    @inbounds productintegrals[ii, ll, p1, p2] = dot(wp1x, pevalx[ll, p2, :]) #p1x .* pll[1][p2].(x)
                end
            end
        end
    end
    return productintegrals
end

"""
    mapfromto(sgfrom, sgto)

Identifies the map for points in grid sgfrom into the grid sgto 
"""
function mapfromto(sgfrom, sgto)
    (fromptids, _, _) = sparsegridpoints(sgfrom.data.terms, sgto.data.terms_to_unique_points_per_dimension, sgto.data.unique_points_per_dimension, sgto.dims)

    terms_to_grid_points = Vector(undef, size(fromptids, 1))
    for ii in 1:size(fromptids, 1)
        terms_to_grid_points[ii] = findfirst(fromptids[ii, :] == uniquerow for uniquerow in eachrow(sgto.data.points_to_unique_indices))
    end
    @debug "Mapped each sparse grid point in sgfrom to sgto"
    return terms_to_grid_points
end

"""
    margin(multi_index_set)
Computes the margin  for a matrix representation of multi-index set.
Each column represents a MI.
"""
function margin(multi_index_set)
    ndim = size(multi_index_set, 1)
    margin = zeros(Int32, ndim, 0)
    for mi in eachcol(multi_index_set)
        for i = 1:ndim
            ei = zeros(Int32, ndim)
            ei[i] = 1
            miplus = mi .+ ei
            if isnothing(findfirst(miplus == mi for mi in eachcol(multi_index_set)))
                margin = hcat(margin, miplus)
            end
        end
    end
    uniqueindices = unique(i -> margin[:, i], 1:size(margin, 2))
    return sortmi(margin[:, uniqueindices])
end

"""
    sortmi(multi_index_set, combination_coeff)

Lexographically sorts a matrix multi index set, and applys the same sorting to combination_coeff.
Each column represents a MI.
"""
function sortmi(multi_index_set, combination_coeff)
    sorted = sortmi([multi_index_set;transpose(combination_coeff)])
    multi_index_set_sorted = sorted[1:end-1,:]
    combination_coeff_sorted = copy(combination_coeff)
    combination_coeff_sorted .= sorted[end,:]
    return (multi_index_set_sorted, combination_coeff_sorted)
end

"""
    sortmi(multi_index_set, combination_coeff)
    
Lexographically sorts the multi index set.
Each column represents a MI.
"""
function sortmi(multi_index_set)
    if ndims(multi_index_set) == 1
        sorted = sort(multi_index_set)
    else
        sorted = sortslices(multi_index_set, dims=2)
    end
    return sorted
end

"""
    reducedmargin(multi_index_set)

Computes the reduced margin for a matrix representation of multi-index set.
Each column represents a MI.
"""
function reducedmargin(multi_index_set)
    mgn = margin(multi_index_set)
    ndim = size(mgn, 1)
    rm = Matrix{Int64}(undef,ndim,0)
    I_ndim = I(ndim)
    for mi in eachcol(mgn)
        in_rm = true

        for ei in eachcol(I_ndim)
            miminus = mi .- ei
            test = in(miminus, eachcol(multi_index_set)) || any(==(0), miminus)
            if test == false
                in_rm = false;
                break
            end
        end
        if in_rm == true
            rm = hcat(rm, mi)
        end
    end
    return sortmi(rm)
end

"""
    createsparsegrid(multi_index_set; rule=Doubling(), knots=CCPoints())

Constructs a SparseGrid for a multi-index set with specified  knots and rules.

# Arguments
- `multi_index_set`: matrix with MI as columns
- `rule`: Level-to-knots rule
- `knots`: Knots type

# Returns
- `sparsegrid`: Constructed sparse grid
"""
function createsparsegrid(multi_index_set; rule=Doubling(), knots=CCPoints())
    multi_index_set = sortmi(multi_index_set)
    maxmi = maximum(multi_index_set)
    dims = size(multi_index_set, 1)

    @debug "Creating sparse grid with parameter dimension "*string(dims)

    if isa(knots, Points)
        knots = fill(knots, dims)
        @debug "Using the same notes in all dimensions"
    end
    if isa(rule, Level)
        rule = fill(rule, dims)
        @debug "Using the same rule in all dimensions"
    end
    try
        @assert (size(rule,1) == size(knots,1))
        @assert (size(rule,1) == dims)
    catch e
        DimensionMismatch("If specified, rule and knots must be the same length as the number of dimensions")
    end
    unique_points_per_dimension = Vector{Vector{Float64}}(undef,dims)
    point_sequences_per_dimension = Vector{Vector{Vector{Float64}}}(undef,dims)
    weight_sequences_per_dimension = Vector{Vector{Vector{Float64}}}(undef,dims)
    terms_to_unique_points_per_dimension = Vector{Dict{Vector{Int},Int}}(undef,dims)
    domain = Vector{Vector{Float64}}(undef,dims)
    for ii = eachindex(knots)
        domain[ii]=knots[ii].domain
        (unique_points_per_dimension[ii], point_sequences_per_dimension[ii], terms_to_unique_points_per_dimension[ii], weight_sequences_per_dimension[ii]) = sparsegridonedimdata(maxmi, knots[ii], rule[ii], domain[ii])
    end
    @debug "Constructed one dimensional grid data"

    (grid, coeff_per_term, combination_coeff, multi_index_set) = sparsegridterms(multi_index_set, point_sequences_per_dimension)
    @debug "Constructed combintion technique data"

    filter_vector = coeff_per_term .!= 0
    grid = grid[filter_vector]
    coeff_per_term = coeff_per_term[filter_vector]

    (points_to_unique_indices, grid_points, terms_to_grid_points) = sparsegridpoints(grid, terms_to_unique_points_per_dimension, unique_points_per_dimension, dims)
    @debug "Created sparse grid term to sparse grid points mappings"

    data = sparse_grid_data{dims}(grid, points_to_unique_indices, terms_to_grid_points, coeff_per_term, terms_to_unique_points_per_dimension, point_sequences_per_dimension, unique_points_per_dimension, weight_sequences_per_dimension)

    sparsegrid = SparseGrid(dims, domain, multi_index_set, combination_coeff, grid_points, knots, rule, data)

    @debug "Created sparse grid"
    return sparsegrid
end

"""
    sparsegridpoints(grid, terms_to_unique_points_per_dimension::Vector{Dict{Vector{Int}, Int}}, unique_points_per_dimension, dims)

Generate sparse grid and associated data.

# Arguments
`grid`: SparseGridTerm vector
`terms_to_unique_points_per_dimension`: Dictionary mapping term to unique 1D grid
`unique_points_per_dimension`: Vector of unique 1d points
`dims`: Dimensionality of sparse grid.

# Returns
- `points_to_unique_indices`: Unique sparse grid in unique 1d point ID representation (i.e. grid represented as N^dims).
- `grid_points`: Unique grid in parameter representation (i.e. grid represneted in parameter domain Gamma)
- `terms_to_grid_points`: Map terms to unique grid indices
"""
function sparsegridpoints(grid, terms_to_unique_points_per_dimension::Vector{Dict{Vector{Int}, Int}}, unique_points_per_dimension, dims)
    # Get the discrete representation of grid as unqiue indices in each dimension
    gridasptsindices = Matrix{Int}(undef, length(grid), dims)
    for ii = eachindex(grid)
        for jj = eachindex(grid[ii].data)
            gridasptsindices[ii, jj] = terms_to_unique_points_per_dimension[jj][grid[ii].data[jj]]
        end
    end

    # Now find unique points
    N = size(gridasptsindices,2)
    row_dict = Dict{NTuple{N, Int}, Int}()
    points_to_unique_indices = NTuple{N, Int}[]
    for i in 1:size(gridasptsindices, 1)
        row = gridasptsindices[i, :]
        tup = ntuple(j -> row[j], dims)
        if !haskey(row_dict, tup)
            row_dict[tup] = length(points_to_unique_indices) + 1
            push!(points_to_unique_indices, tup)
        end
    end

    # # Convert back to matrix
    # Lookup the index for each row in gridasptsindices
    terms_to_grid_points = Vector{Int}(undef, size(gridasptsindices,1))
    for ii in 1:size(gridasptsindices,1)
        row = gridasptsindices[ii,:]
        tup = ntuple(j -> row[j], dims)
        terms_to_grid_points[ii] = row_dict[tup]
    end

    grid_points = Matrix{Float64}(undef,length(points_to_unique_indices),dims)
    for ii = eachindex(points_to_unique_indices)
        for jj = 1:dims
            grid_points[ii, jj] = unique_points_per_dimension[jj][points_to_unique_indices[ii][jj]]
        end
    end

    return (points_to_unique_indices, grid_points, terms_to_grid_points)
end

"""
    sparsegridterms(multi_index_set, point_sequences_per_dimension; combination_coeff=nothing)

Builds the sparse grid terms in level and point indices, and computes corresponding combination technique coeffs.
"""
function sparsegridterms(multi_index_set, point_sequences_per_dimension; combination_coeff=nothing)
    if isnothing(combination_coeff)
        combination_coeff = computesparsegridc(multi_index_set)
    end
    n = size(multi_index_set,1)
    nterms_batch = 1000
    terms = Vector{SparseTerm{n}}(undef,nterms_batch)
    coeff_per_term = Vector{Int}(undef,nterms_batch)

    count = 1
    length_terms = nterms_batch

    for (miind, mi) in enumerate(eachcol(multi_index_set))
        if true #combination_coeff[miind] != 0
            indexarray = Vector{Vector{Tuple{Int, Int}}}(undef, 0)  # outer vector holds vectors of tuples
            for (ind,i) in enumerate(mi)
                indexarrayi = Vector{Tuple{Int, Int}}()  # vector of tuples, no fixed size needed here
                if isa(point_sequences_per_dimension[1][1], Number)
                    for j in eachindex(point_sequences_per_dimension[i])
                        push!(indexarrayi, (i, j))   # use tuple, immutable and stack-allocated
                    end
                else
                    for j in eachindex(point_sequences_per_dimension[ind][i])
                        push!(indexarrayi, (i, j))   # tuple again
                    end
                end
                push!(indexarray, indexarrayi)
            end

            # for row in collect(Iterators.product((indexarray[i] for i in eachindex(indexarray))...))
            for row in Iterators.product((indexarray[i] for i in eachindex(indexarray))...)
                t = build_sparse_term(Val(n), row) 
                terms[count] = t
                # push!(coeff_per_term, combination_coeff[miind])
                coeff_per_term[count] = combination_coeff[miind]
                count +=1
                if count > length_terms
                    resize!(terms, length_terms + nterms_batch)
                    resize!(coeff_per_term, length(coeff_per_term) + nterms_batch)
                    length_terms += nterms_batch
                end
            end
        end
    end
    # Now trim terms and coeff_per_count
    resize!(terms, count - 1)
    resize!(coeff_per_term, count - 1)
    return (terms, coeff_per_term, combination_coeff, multi_index_set)
end

"""
    build_sparse_term(::Val{n}, row)

Constructor for SparseTerm.
    Optimised code.
"""
@generated function build_sparse_term(::Val{n}, row) where n
    # Unroll tuple construction for performance
    return :(SparseTerm{$n}(tuple($(map(i -> :(SVector{2,Int}(row[$i][1], row[$i][2])), 1:n)...))))
end

"""
   sparsegridonedimdata(maxmi, knots, rule, domain)

For set of knots and level-to-knots rule, builds point and weight sequences.
"""
@noinline function sparsegridonedimdata(maxmi, knots, rule, domain)
    unique_points_per_dimension = Vector{Float64}(undef, 0)
    point_sequences_per_dimension = Vector{Vector{Float64}}(undef, maxmi)
    weight_sequences_per_dimension = Vector{Vector{Float64}}(undef, maxmi)
    terms_to_unique_points_per_dimension = Dict{Vector{Int},Int}()

    for ii = 1:maxmi
        # Detect special case of fidelity
        if isa(knots, FidelityPoints)
            # This must be treated differently, with no lev-2-knots rule.
            point_sequences_per_dimension[ii], weight_sequences_per_dimension[ii] = knots(ii)
        else
            point_sequences_per_dimension[ii], weight_sequences_per_dimension[ii] = knots(rule(ii))
        end

        # Build points map
        for (jj, p) in enumerate(point_sequences_per_dimension[ii])
            if length(unique_points_per_dimension) == 0
                ind = nothing
            else
                ind = findfirst(x -> isapprox(p, x), unique_points_per_dimension)
            end

            if !isnothing(ind)
                terms_to_unique_points_per_dimension[[ii; jj]] = ind
            else
                push!(unique_points_per_dimension, p)
                terms_to_unique_points_per_dimension[[ii; jj]] = length(unique_points_per_dimension)
            end
        end
    end
    return (unique_points_per_dimension, point_sequences_per_dimension, terms_to_unique_points_per_dimension, weight_sequences_per_dimension)
end

"""
    interpolateonsparsegrid(sparsegrid, fongrid, targetpoints; evaltype=nothing)

Interpolates the evaluations fongrid with sparsegrid onto targetpoints

# Arguments
- `sparsegrid`: SparseGrid structure
- `fongrid`: Vector of function evaluations at get_grid_points(sparsegrid)
- `targetpoints`: Vector of targetpoints

# Returns
- `feval`: Vector of interpolation to targetpoints
"""
function interpolateonsparsegrid(sparsegrid, fongrid, targetpoints; evaltype=nothing) 
    # Set up return vector
    if isempty(fongrid)
        zero_elem = zero(evaltype)
    else
        zero_elem = zero(fongrid[1])
    end
    feval = Vector{typeof(zero_elem)}(undef, length(targetpoints))
    copyable = try 
        copy(zero_elem) 
        true
     catch 
        false 
    end
    if copyable
        for i in eachindex(feval)
            feval[i] = copy(zero_elem)
        end
    else
        # Assume produces a new object
        @warn "Object cannot be copied, assuming scalar multiplication creates new object"
        for i in eachindex(feval)
            feval[i] = 0.0*zero_elem
        end
    end
    terms = sparsegrid.data.terms
    coeff_per_term = sparsegrid.data.coeff_per_term
    maprowtouniquept = sparsegrid.data.terms_to_grid_points
    nparams = sparsegrid.dims

    ntarget = length(targetpoints)
    targetpoints_matrix = copy(reduce(hcat, targetpoints)')
    acc = Vector{Float64}(undef,ntarget)
    f_placeholder = Vector{Float64}(undef,ntarget)
    # cweightedpolyprod = zeros(length(feval),size(terms,1))
    cweightedpolyprod_i = zeros(Float64,length(feval),1)

    # Runtime check if `fongrid` supports indexing and mutation
    process_elementwise = try
        # Test if fongrid supports indexing
        _ = fongrid[1]
        # Test if feval[1] supports eachindex
        eachindex(fongrid[1][1])
        true
    catch
        @warn "Object cannot be indexed, this may introduce intermediate allocation"
        false
    end

    # Precompute lagrange interpolation once to unique points
    dims = sparsegrid.dims
    maxmi = maximum(sparsegrid.multi_index_set,dims=2)
    lagrange_evaluation_precompute = Vector{Vector{Vector{Vector{Float64}}}}(undef,dims)
    for j in 1:dims
        targetpoints_j = @view targetpoints_matrix[:,j]
        maxmi_j = maxmi[j]
        lagrange_evaluation_precompute[j] = Vector{Vector{Vector{Float64}}}(undef,maxmi_j)
        for r1 = 1:maxmi_j
            pts = sparsegrid.data.point_sequences_per_dimension[j][r1]
            n_pts = length(pts)
            lagrange_evaluation_precompute[j][r1] = Vector{Vector{Float64}}(undef,n_pts)
            for r2 = 1:n_pts
                lagrange_evaluation_precompute[j][r1][r2] = Vector{Float64}(undef,length(targetpoints_j))
                lagrange_evaluation!(lagrange_evaluation_precompute[j][r1][r2],
                pts,r2,targetpoints_j)
            end
        end
    end

    # Iterate through terms updating the evaluataions at target points
    for i in eachindex(terms)
        row = terms[i].data
        @inbounds begin
            acc .= 1.0
            # Extract appropriate evaluations and accumulate
            for j in eachindex(row)
                @inbounds targetpoints_j = @view targetpoints_matrix[:,j]
                r1 = row[j][1]
                r2 = row[j][2]
                f_placeholder = lagrange_evaluation_precompute[j][r1][r2]
                acc .*= f_placeholder
            end
            mul!(cweightedpolyprod_i,coeff_per_term[i],acc)
        end
        # Update evaluations with term
        fused_update!(feval, fongrid[maprowtouniquept[i]], cweightedpolyprod_i,process_elementwise)
    end
    return feval
end

"""
    fused_update!(feval,fsrc,cweightedpolyprod_i,process_elementwise)

Update feval based upon a function evaluation and appropriate scaling for term and targetpoints
This aims to avoid allocation if possible (when process_elementwise==true)
"""
function fused_update!(feval,fsrc,cweightedpolyprod_i,process_elementwise)
    scale = cweightedpolyprod_i
    @inbounds for kk in eachindex(feval, scale)
        y = feval[kk]
        α = scale[kk]

        if process_elementwise
            # Non-allocating elementwise update
            @inbounds @simd for j in eachindex(y, fsrc)
                y[j] += α * fsrc[j]
            end
        else
            # Fallback to allocating
            feval[kk] = y + α * fsrc
        end
    end
end

"""
    integrateonsparsegrid(sparsegrid, fongrid, precompute; evaltype=nothing)

NOTE: This has been replaced by precomputed quadrature weights.
Computes integral of approximation defined via sparsegrid and fongrid.

Requires precomputed integrals of pairwise products of interpolation polynomials.
"""
function integrateonsparsegrid(sparsegrid, fongrid, precompute; evaltype=nothing)
    if isempty(fongrid)
        if isnothing(evaltype)
            error("Must input evaltype if using zero operator")
        end
        fongrid[1] = 0.0*zero(evaltype)
    else
        evaltype = typeof(fongrid[1])
    end
    terms = sparsegrid.data.terms
    coeff_per_term = sparsegrid.data.coeff_per_term
    maprowtouniquept = sparsegrid.data.terms_to_grid_points
    nterms = size(terms,1)
    nparams = sparsegrid.dims
    if length(precompute.productintegrals) == 1
        productintegrals = Vector{Any}(undef,nparams) 
        for ii = 1:nparams
            productintegrals[ii] = precompute.productintegrals[1]
        end
    else
        productintegrals = precompute.productintegrals 
    end

    val = Vector{Float64}(undef,nterms)
    for (i, row) = enumerate(eachrow(terms))
        val[i] = coeff_per_term[i] * prod(productintegrals[j][row[1][j][1], 1, row[1][j][2], 1] for j in eachindex(row[1]))
    end
    @debug "Computed integral on sparse grid"
    feval = sum(val[i] * fongrid[maprowtouniquept[i]] for i=1:nterms)
    return feval
end

"""
    computepairwisenorms(fongrid; product=nothing)

Computes the product(x,y) for pairs of function evaluations.
Default product is x.*y.

Can be used for example with a mass matrix product = (x,y) -> sqrt(x' * Q * y)
"""
function computepairwisenorms(fongrid; product=nothing)
    if isnothing(product)
        product = (x,y) -> x.*y
    end
    ngrid = length(fongrid)
    pairwisenorms = Matrix{typeof(product(fongrid[1], fongrid[1]))}(undef, (ngrid, ngrid))
    for i1 = 1:ngrid, i2 = 1:ngrid
            pairwisenorms[i1, i2] = product(fongrid[i1], fongrid[i2])
    end
    return pairwisenorms
end

"""
    L2onsparsegrid(sparsegrid, fongrid, precompute; product=nothing, pairwisenorms=nothing, returnpairwisenorms=false)

Computes weighted L^2 integral based upon precomputed pairwise products of interpolation polynomials and pairwise products.
"""
function L2onsparsegrid(sparsegrid, fongrid, precompute; product=nothing, pairwisenorms=nothing, returnpairwisenorms=false)
    if isnothing(pairwisenorms)
        @debug "Must compute pairwise norms - may be possible to precompute for efficiency"
        pairwisenorms = computepairwisenorms(fongrid, product=product)
    end

    terms = sparsegrid.data.terms
    coeff_per_term = sparsegrid.data.coeff_per_term
    maprowtouniquept = sparsegrid.data.terms_to_grid_points
    nparams = sparsegrid.dims
    if length(precompute.productintegrals) == 1
        productintegrals = Vector{Any}(undef,nparams) 
        for ii = 1:nparams
            productintegrals[ii] = precompute.productintegrals[1]
        end
    else
        productintegrals = precompute.productintegrals 
    end
    feval = zeros(ComplexF64, size(pairwisenorms[1,1]))

    nonzero_rows = [any(!iszero(pairwisenorms[maprowtouniquept[i], :])) for i in 1:size(terms,1)]
    valstomultiply = ones(size(terms,1),size(terms,1))
    n = length(terms[1].data)

    a = Vector{Int64}(undef,2)
    b = Vector{Int64}(undef,2)

    for (i, rowi) = enumerate(terms), (j, rowj) = enumerate(terms)
        @inbounds if nonzero_rows[i] == false
            continue
        end
        @inbounds if nonzero_rows[j] == false
            continue
        end
        rowi1 = rowi.data
        rowj1 = rowj.data
        @inbounds for k in eachindex(rowi1)
            a .= rowi1[k]
            b .= rowj1[k]
            @fastmath valstomultiply[i,j] *= productintegrals[k][a[1],b[1],a[2],b[2]]
        end
        @inbounds feval .= feval .+ coeff_per_term[i] * coeff_per_term[j] * valstomultiply[i,j] * pairwisenorms[maprowtouniquept[i], maprowtouniquept[j]]
    end

    feval = sqrt.(abs.(feval)) # abs prevents very small negative numbers (≈ 0) giving error

    if returnpairwisenorms == true
        return (feval, pairwisenorms)
    else
        return feval
    end
end

"""
    computesparsegridc_old(I)

Computes combination technique coefficients for multi-index matrix I
This is a slower algorithm but does not depend on admissibility of I
"""
function computesparsegridc_old(I)
    jjs = generatebinaryvectors(size(I, 1))
    combination_coeff = Vector{Integer}(undef, size(I, 2))
    combination_coeff .= 0
    iiplusjjs = similar(jjs)
    for (ii, i) = enumerate(eachcol(I))
        iiplusjjs = jjs .+ i
        for (jj, ij) = enumerate(eachcol(iiplusjjs))
            if any(ij == combination_coeff for combination_coeff in eachcol(I))
                combination_coeff[ii] = combination_coeff[ii] .+ (-1)^norm(jjs[:, jj], 1)
            end
        end
    end
    return combination_coeff
end

"""
    computesparsegridc(I)

Computes combination technique coefficients for multi-index matrix I.
Branches depending on admissibility.
"""
function computesparsegridc(I)
    (adm, mi_missing) = check_admissibility(MISet([vec(collect(c)) for c in eachcol(I)]))
    if adm
        coeff = computesparsegridc_admissible(I)
    else
        # Slower
        coeff = computesparsegridc_old(I)
    end
end

"""
    computesparsegridc(I)

Computes combination technique coefficients for an admissible multi-index matrix I.

Based upon the Sparse Grids Matlab Kit code.
https://sites.google.com/view/sparse-grids-kit
@article{piazzola.tamellini:SGK,
 author = {Piazzola, C. and Tamellini, L.},
 title  = {{Algorithm 1040: The Sparse Grids Matlab Kit - a Matlab implementation of sparse grids for high-dimensional function approximation and uncertainty quantification}},
 journal= {ACM Transactions on Mathematical Software},
 year   = {2024},
 volume = {50},
 number = {1},
 doi = {10.1145/3630023}
}
"""
function computesparsegridc_admissible(I)
    # Based upon SGMK
    nn = size(I, 2)
    coeff = ones(Int, nn)

    bookmarks = unique(i -> I[1,i], 1:size(I,2))
    bk = [bookmarks[3:end]'.-1 nn nn];
    for i in 1:nn
        cc = I[:, i]
        range = bk[cc[1]]
        for j in (i+1):range
            d = I[:, j] .- cc
            if maximum(d) <= 1 && minimum(d) >= 0
                coeff[i] += (-1)^sum(d)
            end
        end
    end
    return coeff
end

"""
    generatebinaryvectors(n)
    Generates all vectors in {0,1}^n
"""
function generatebinaryvectors(n)
    vectors = Matrix{Integer}(undef, n, 0)
    for i in 0:2^n-1
        binary_rep = reverse(digits(i, base=2, pad=n))
        vectors = hcat(vectors, binary_rep)
    end
    return vectors
end

"""
    createsmolyakmiset(dim::Int, level::Int)

Creates Smolyak-type multi-index set in dim dimensions.
That is, || α ||_1 <= dim + level
"""
function createsmolyakmiset(dim::Int, level::Int)
    multi_index_list = Vector{Vector{Int}}()
    current_index = Vector{Int}(undef, dim)

    total_order_min = dim
    total_order_max = dim + level

    # Loop through permissible || α ||_1
    for total_order in total_order_min:total_order_max
        generate_compositions!(multi_index_list, current_index, total_order, dim, 1)
    end

    return sortmi(hcat(multi_index_list...))
end

"""
    generate_compositions!(
    results::Vector{Vector{Int}},
    current::Vector{Int},
    remaining::Int,
    parts::Int,
    idx::Int)

Generates compositions of remaining in the indices idx to end

Used to recursively generate the Smolyak MI set
"""
function generate_compositions!(
    results::Vector{Vector{Int}},
    current::Vector{Int},
    remaining::Int,
    parts::Int,
    idx::Int,
)
    # If we are at the last slot, fill it with the remainder and record the composition
    if idx == parts
        current[idx] = remaining
        push!(results, copy(current))
        return
    end
    # We must leave at least 1 for each remaining slot
    limit = remaining - (parts - idx)
    for val in 1:limit
        current[idx] = val
        generate_compositions!(
            results,
            current,
            remaining - val,
            parts,
            idx + 1,
        )
    end
end

"""
    createtensormiset(dim::Int, level::Int)

Creates tensor product multi-index set
"""
function createtensormiset(dim::Int, level::Int)
    total_points = (level + 1)^dim
    output = Matrix{Int}(undef, dim, total_points)
    current_index = fill(1, dim)

    for col in 1:total_points
        @inbounds for d in 1:dim
            output[d, col] = current_index[d]
        end

        # Increment multi-index like a counter
        for d in dim:-1:1
            current_index[d] += 1
            if current_index[d] <= level + 1
                break
            else
                current_index[d] = 1
            end
        end
    end

    return sortmi(output)
end

"""
    is_leq(α::Vector, β::Vector)
Componentwise leq for vectors
"""
function is_leq(α::Vector, β::Vector)
    # Check if β ≤ α element-wise
    return all(β[i] ≤ α[i] for i in 1:length(α))
end

"""
    downwards_closed_set(input_matrix::Matrix{Int})

Tests admissibility for matrix format multi-index set
"""
function downwards_closed_set(input_matrix::Matrix{Int})
    num_indices, num_columns = size(input_matrix)

    # Use a Set to store seen multi-indices as NTuples (cheaper lookup)
    seen = Set{NTuple{typeof(input_matrix[1,1]), Int}}()
    result = NTuple{num_indices, Int}[]

    for col_idx in 1:num_columns
        α = view(input_matrix, :, col_idx)
        # Generate all β ≤ α
        for β_tuple in Iterators.product((0:α[i] for i in 1:num_indices)...)
            if !in(β_tuple, seen)
                push!(result, β_tuple)
                push!(seen, β_tuple)
            end
        end
    end

    # Convert result to a matrix
    output_matrix = hcat(result...)
    return sortmi(output_matrix)
end

"""
    check_unique_pts(pts)
Ensures pts is a unique sequence
"""
function check_unique_pts(pts)
    @inbounds for i in 1:length(pts)-1
        pi = pts[i]
        @inbounds for j in i+1:length(pts)
            if pi == pts[j]
                error("Duplicate points detected: pts[$i] == pts[$j] == $(pts[j])")
            end
        end
    end
    return nothing
end

"""
    lagrange_evaluation!(y, pts, ii::Int, targetpoints)

Evaluates Lagrange polynomial defined by pts for pts(ii) at targetpoints
"""
@inline function lagrange_evaluation!(y, pts, ii::Int, targetpoints)
    @inbounds begin
        denom = 1.0
        x_ii = pts[ii]
        for j in eachindex(pts)
            j == ii && continue
            denom *= x_ii - pts[j]
        end

        pts_set = Set(pts)

        for k in eachindex(targetpoints)
            x = targetpoints[k]

            # Check if x matches any interpolation node exactly
            exact_match = x ∈ pts_set
            if exact_match
                for j in eachindex(pts)
                    if x == pts[j]
                        if j == ii
                            y[k] = 1.0
                        else
                            y[k] = 0.0
                        end
                        break
                    end
                end
                continue
            end

            # Compute numerator product normally
            numer = 1.0
            for j in eachindex(pts)
                j == ii && continue
                numer *= x - pts[j]
            end
            y[k] = numer / denom
        end
    end
end