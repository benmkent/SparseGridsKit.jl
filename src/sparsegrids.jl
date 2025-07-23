using LinearAlgebra
using StaticArrays
using Polynomials

import Base.:+
import Base.:-

mutable struct sparse_grid_data
    terms::Vector{Any} # Terms in the sparse grid
    points_to_unique_indices::Matrix{Int} # Mapping from terms to unique points
    terms_to_grid_points::Vector{Int} # Mapping from terms to grid points
    coeff_per_term::Vector{Int} # Combination coefficient per term
    terms_to_unique_points_per_dimension::Vector{Dict{Vector{Int},Int}} # Mapping from terms to unique points per dimension
    point_sequences_per_dimension::Vector{Vector{Vector{Float64}}}
    unique_points_per_dimension::Vector{Vector{Float64}} # Unique points per dimension
end

mutable struct SparseGrid
    dims::Int # Number of dimensions
    domain::Vector{Vector{Real}} # Grid domain
    multi_index_set::Matrix{Int} # Matrix of multi-index columns [β_1 β_2 ...]
    combination_coeff::Vector{Int} # Combination technique coefficients
    grid_points::Matrix{Real} # Sparse grid
    quadrature_weights::Vector{Float64} # Quadrature weights
    knots::Vector{Any}
    rules::Vector{Any}
    data::sparse_grid_data # Structure for underlying construction data
end

function SparseGrid(dims, domain, multi_index_set, combination_coeff, grid_points, knots, rules, data)
    sg = SparseGrid(dims, domain, multi_index_set, combination_coeff, grid_points, [], knots, rules, data)
    return sg
end

mutable struct SparseGridApproximation
    sg::SparseGrid
    fongrid::Vector
    domain::Vector
end

function SparseGridApproximation(sg, fongrid)
    SparseGridApproximation(sg,fongrid,sg.domain)
end

function (sga::SparseGridApproximation)(x)
    @assert length(x) == sga.sg.dims
    return interpolate_on_sparsegrid(sga.sg, sga.fongrid, [x])[1]
end

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

    # Clean zero c points
    #zeroindices = sg.combination_coeff .!= 0
    #combination_coeff = sg.combination_coeff[zeroindices]
    #multi_index_set = sg.multi_index_set[:, zeroindices]

    (multi_index_set,combination_coeff) = sortmi(sg.multi_index_set, sg.combination_coeff)

    (terms, coeff_per_term, combination_coeff, multi_index_set) = sparsegridterms(multi_index_set, sg.data.point_sequences_per_dimension, combination_coeff=sg.combination_coeff)

    (points_to_unique_indices, grid_points, terms_to_grid_points) = sparsegridpoints(terms, sg.data.terms_to_unique_points_per_dimension, sg.data.unique_points_per_dimension, size(sg.multi_index_set, 1))

    data = sparse_grid_data(terms, points_to_unique_indices, terms_to_grid_points, coeff_per_term, terms_to_unique_points_per_dimension, sg.data.point_sequences_per_dimension, sg.data.unique_points_per_dimension)

    sparsegrid = SparseGrid(dims, sg.domain, multi_index_set, combination_coeff, grid_points, sg.knots, sg.rules, data)
    return sparsegrid
end

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
        (unique_points_per_dimension, point_sequences_per_dimension, terms_to_unique_points_per_dimension) = sparsegridonedimdata(maxmi, knots[ii], rule[ii], knots[ii].domain)
        productintegrals[ii] = sparsegridproductintegrals(point_sequences_per_dimension, maxmi, knots[ii], rule[ii])
    end
    @debug "Computed weighted L^2 product integrals in each dimension"
    return (; productintegrals)
end

function sparsegridproductintegrals(point_sequences_per_dimension, maxmi, knots, rule)
    @debug "Computing product integrals"

    # Check uniqueness of points
    for pts in point_sequences_per_dimension
        check_unique_pts(pts)
    end

    # Allocate arrays
    productintegrals = zeros(Float64, maxmi, maxmi, length(point_sequences_per_dimension[end]), length(point_sequences_per_dimension[end]))

    #x, w = gausslegendre(2 * rule(maxmi))
    x,w = knots(2*rule(maxmi))
    wp1x = similar(x)

    pevalx = zeros(Float64, maxmi, length(point_sequences_per_dimension[end]), length(x))
    for ii = 1:maxmi
        # for jj = 1:length(point_sequences_per_dimension[ii])
        for jj = 1:length(point_sequences_per_dimension[ii])
            # pevalx[ii, jj, :] .= polyperlevel[ii][jj].(x)
            write_location = @view pevalx[ii, jj, :]
            pts_input = point_sequences_per_dimension[ii]
            lagrange_evaluation!(write_location, pts_input, jj, x)
        end
    end

    #productintegralsmatrix = Matrix{Float64}(undef, length(polyperlevel[end]), length(polyperlevel[end]))
    for ii = 1:maxmi
        # Compute integrals
        #productintegralslevel = []
        #pii = view(polyperlevel,ii);
        for ll = 1:maxmi
            #pll = view(polyperlevel,ll);s
            @inbounds for p1 in eachindex(point_sequences_per_dimension[ii])
                #p1x .= pii[1][p1].(x)
                #p1x .= w.* p1x
                wp1x .= w .* pevalx[ii, p1, :]
                @inbounds for p2 in eachindex(point_sequences_per_dimension[ll])
                    #productintegralsmatrix[p1, p2] = sum(w .* polyperlevel[ii][p1].(x) .* polyperlevel[ll][p2].(x))
                    @inbounds productintegrals[ii, ll, p1, p2] = dot(wp1x, pevalx[ll, p2, :]) #p1x .* pll[1][p2].(x)
                    #display([ii,ll,p1,p2,sum(w .* polyperlevel[ii][p1].(x) .* polyperlevel[ll][p2].(x))])
                end
            end
            #push!(productintegralslevel, productintegralsmatrix[1:length(polyperlevel[ii]), 1:length(polyperlevel[ll])])
        end
        #push!(productintegrals, productintegralslevel)
    end
    return productintegrals
end

function mapfromto(sgfrom, sgto)
    (fromptids, ~, ~) = sparsegridpoints(sgfrom.data.terms, sgto.data.terms_to_unique_points_per_dimension, sgto.data.unique_points_per_dimension, sgto.dims)

    terms_to_grid_points = Vector(undef, size(fromptids, 1))
    for ii in 1:size(fromptids, 1)
        terms_to_grid_points[ii] = findfirst(fromptids[ii, :] == uniquerow for uniquerow in eachrow(sgto.data.points_to_unique_indices))
    end
    @debug "Mapped each sparse grid point in sgfrom to sgto"
    return terms_to_grid_points
end


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

function sortmi(multi_index_set, combination_coeff)

    sorted = sortmi([multi_index_set;transpose(combination_coeff)])
    multi_index_set_sorted = sorted[1:end-1,:]
    combination_coeff_sorted = copy(combination_coeff)
    combination_coeff_sorted .= sorted[end,:]
    return (multi_index_set_sorted, combination_coeff_sorted)
end

function sortmi(multi_index_set)
    if ndims(multi_index_set) == 1
        sorted = sort(multi_index_set)
    else
        sorted = sortslices(multi_index_set, dims=2)
    end
    return sorted
end

function reducedmargin(multi_index_set)
    mgn = margin(multi_index_set)
    ndim = size(mgn, 1)
    rm = Matrix{Int64}(undef,ndim,0)
    I_ndim = I(ndim)
    for mi in eachcol(mgn)
        in_rm = true

        for ei in eachcol(I_ndim)
            miminus = mi .- ei
            
            # test = any(v -> v == miminus, eachcol(multi_index_set)) || any(miminus .== 0)
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
    terms_to_unique_points_per_dimension = Vector{Dict{Vector{Int},Int}}(undef,dims)
    domain = Vector{Vector{Float64}}(undef,dims)
    for ii = eachindex(knots)
        domain[ii]=knots[ii].domain
        (unique_points_per_dimension[ii], point_sequences_per_dimension[ii], terms_to_unique_points_per_dimension[ii]) = sparsegridonedimdata(maxmi, knots[ii], rule[ii], domain[ii])
    end
    @debug "Constructed one dimensional grid data"

    (grid, coeff_per_term, combination_coeff, multi_index_set) = sparsegridterms(multi_index_set, point_sequences_per_dimension)
    @debug "Constructed combintion technique data"

    filter_vector = coeff_per_term .!= 0
    grid = grid[filter_vector]
    coeff_per_term = coeff_per_term[filter_vector]

    (points_to_unique_indices, grid_points, terms_to_grid_points) = sparsegridpoints(grid, terms_to_unique_points_per_dimension, unique_points_per_dimension, dims)
    @debug "Created sparse grid term to sparse grid points mappings"

    data = sparse_grid_data(grid, points_to_unique_indices, terms_to_grid_points, coeff_per_term, terms_to_unique_points_per_dimension, point_sequences_per_dimension, unique_points_per_dimension)

    sparsegrid = SparseGrid(dims, domain, multi_index_set, combination_coeff, grid_points, knots, rule, data)

    @debug "Created sparse grid"
    return sparsegrid
end

function sparsegridpoints(grid, terms_to_unique_points_per_dimension, unique_points_per_dimension, dims)
    if ~isa(terms_to_unique_points_per_dimension,Vector)
        commonknots = true
    else
        commonknots = false
    end

    gridasptsindices = Matrix{Integer}(undef, length(grid), dims)
    for ii = eachindex(grid)
        for jj = eachindex(grid[ii])
            if commonknots == true
                gridasptsindices[ii, jj] = terms_to_unique_points_per_dimension[grid[ii][jj]]
            else
                gridasptsindices[ii, jj] = terms_to_unique_points_per_dimension[jj][grid[ii][jj]]
            end
        end
    end

    N = size(gridasptsindices,2)

    # Find unique points
    uniquerows = unique(i -> gridasptsindices[i, :], 1:size(gridasptsindices, 1))
    points_to_unique_indices = gridasptsindices[uniquerows, :]

    # terms_to_grid_points = Vector{Integer}(undef, size(gridasptsindices, 1))
    # for ii in 1:size(gridasptsindices, 1)
    #     # terms_to_grid_points[ii] = findfirst(gridasptsindices[ii, :] == uniquerow for uniquerow in eachrow(points_to_unique_indices))
    #     terms_to_grid_points[ii] = findfirst(row -> gridasptsindices[ii, :] == row, eachrow(points_to_unique_indices))
    # end


    # Build dictionary mapping each unique row (as Tuple) to its index
    row_dict = Dict{NTuple{N, Int}, Int}()

    for (idx, row) in enumerate(eachrow(points_to_unique_indices))
        row_dict[Tuple(row)] = idx
    end

    # Lookup the index for each row in gridasptsindices
    terms_to_grid_points = Vector{Int}(undef, size(gridasptsindices, 1))
    for ii in 1:size(gridasptsindices, 1)
        terms_to_grid_points[ii] = row_dict[Tuple(gridasptsindices[ii, :])]
    end

    grid_points = similar(points_to_unique_indices, Float64)
    for ii = 1:size(points_to_unique_indices, 1)
        for jj = 1:size(points_to_unique_indices, 2)
            if commonknots == true
                grid_points[ii, jj] = unique_points_per_dimension[points_to_unique_indices[ii, jj]]
            else
                grid_points[ii, jj] = unique_points_per_dimension[jj][points_to_unique_indices[ii, jj]]
            end
        end
    end
    return (points_to_unique_indices, grid_points, terms_to_grid_points)
end

function sparsegridterms(multi_index_set, point_sequences_per_dimension; combination_coeff=nothing)
    if isnothing(combination_coeff)
        combination_coeff = computesparsegridc(multi_index_set)
    end
    n = size(multi_index_set,1)
    terms = []
    coeff_per_term = []
    for (miind, mi) in enumerate(eachcol(multi_index_set))
        if true #combination_coeff[miind] != 0
            indexarray = Vector{Vector{Vector{Int}}}(undef, 0)
            for (ind,i) in enumerate(mi)
                indexarrayi = Vector{Vector{Int}}(undef, 0)
                if isa(point_sequences_per_dimension[1][1],Number)
                    for j = eachindex(point_sequences_per_dimension[i])
                        push!(indexarrayi, [i, j])
                    end
                else
                    for j = eachindex(point_sequences_per_dimension[ind][i])
                        push!(indexarrayi, [i, j])
                    end   
                end
                push!(indexarray, collect(indexarrayi))
            end

            for row in collect(Iterators.product((indexarray[i] for i in eachindex(indexarray))...))
                push!(terms,SVector{n}(row))
                push!(coeff_per_term, combination_coeff[miind])
            end
        end
    end
    # Clean up multi_index_set and combination_coeff
    #zeroindices = combination_coeff .!= 0
    #multi_index_set = multi_index_set[:, zeroindices]
    #combination_coeff = combination_coeff[zeroindices]
   

    return (terms, coeff_per_term, combination_coeff, multi_index_set)
end

@noinline function sparsegridonedimdata(maxmi, knots, rule, domain)
    unique_points_per_dimension = Vector{Real}(undef, 0)
    point_sequences_per_dimension = Vector{Vector{Real}}(undef, maxmi)
    #poly = Vector{Vector{Polynomial{Float64,:x}}}(undef, maxmi)
    terms_to_unique_points_per_dimension = Dict{Vector{Int},Int}()

    # Based on knots define appropriate ApproxFun space for underlying polynomials
    # space = spaceselector(domain)

    for ii = 1:maxmi
        # Detect special case of fidelity
        if isa(knots, FidelityPoints)
            point_sequences_per_dimension[ii], w = knots(ii)
        else
            point_sequences_per_dimension[ii], w = knots(rule(ii))
        end

        # Build points map
        for (jj, p) in enumerate(point_sequences_per_dimension[ii])
            if length(unique_points_per_dimension) == 0
                ind = nothing
            else
                # ind = findfirst(isapprox(p, unique_points_per_dimension[i]) for i in eachindex(unique_points_per_dimension))
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
    return (unique_points_per_dimension, point_sequences_per_dimension, terms_to_unique_points_per_dimension)
end

function interpolateonsparsegrid(sparsegrid, fongrid, targetpoints; evaltype=nothing)  
    if isempty(fongrid)
        if isnothing(evaltype)
            error("Must input evaltype if using zero operator")
        end
        feval = [0.0*fongrid[1] for _ in 1:length(targetpoints)]
        return feval
    else
        evaltype = typeof(fongrid[1])
        feval = [0.0*fongrid[1] for _ in 1:length(targetpoints)]
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
    cweightedpolyprod_i = zeros(length(feval),1)
    temp_mul = Ref(copy(fongrid[1]))
    temp_mul[] = 0.0* fongrid[1]
    using_floats = typeof(fongrid[1]) <: Number
    #for k = 1:length(feval)
        for i in eachindex(terms)
        #targetpt = targetpoints[k]
        #for (i, row) = enumerate(eachrow(terms))
            #val = c[i]*prod([polyperlevel[row[1][j][1]][row[1][j][2]](newpt[1][j]) for j in eachindex(row[1])])
            # @inbounds  cweightedpolyprod[k,i] = coeff_per_term[i] * prod(poly[j][row[1][j][1]][row[1][j][2]](targetpt[j]) for j in eachindex(row[1]))
            row = @view terms[i,:]
            @inbounds begin
                acc .= 1.0
                for j in eachindex(row[1])
                    # f = poly[j][row[1][j][1]][row[1][j][2]]
                    # f_placeholder = f.([targetpoints[k][j] for k = 1:ntarget])
                    # targetpoints_j = [targetpoints[k][j] for k = 1:ntarget]
                    @inbounds targetpoints_j = @view targetpoints_matrix[:,j]
                    lagrange_evaluation!(f_placeholder,sparsegrid.data.point_sequences_per_dimension[j][row[1][j][1]],row[1][j][2], targetpoints_j)
                    acc .= acc .* f_placeholder
                end
                # cweightedpolyprod[:,i] = coeff_per_term[i] * acc
                cweightedpolyprod_i .= coeff_per_term[i] * acc
            end
            # feval += cweightedpolyprod[:,i] * fongrid[maprowtouniquept[i]]
            # if using_floats
                # for kk = eachindex(feval)
                #     temp_mul[] = fongrid[maprowtouniquept[i]]
                #     temp_mul[] = temp_mul[] * cweightedpolyprod_i[kk]
                #     @inbounds feval[kk] = feval[kk] + temp_mul[]
                # end
                for kk = eachindex(feval)
                    @inbounds copyto!(temp_mul[], fongrid[maprowtouniquept[i]])
                    @inbounds for j in eachindex(temp_mul[])
                        temp_mul[][j] *= cweightedpolyprod_i[kk]
                    end
                    @inbounds feval[kk] .+= temp_mul[]
                end
            # else
            #     for kk = eachindex(feval)
            #         temp_mul[] = fongrid[maprowtouniquept[i]]
            #         temp_mul[] *= cweightedpolyprod_i[kk]
            #         @inbounds feval[kk] += temp_mul[]
            #     end
            # end
                # feval += cweightedpolyprod_i .* fongrid[maprowtouniquept[i]]
        end
    #end
    # @inbounds @fastmath feval = [sum(cweightedpolyprod[k,i] * fongrid[maprowtouniquept[i]] for i in 1:size(terms,1)) for k=1:length(feval)]

    # mul!(feval, cweightedpolyprod, fongrid[maprowtouniquept])

    return feval
end

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

function computepairwisenorms(fongrid; product=nothing)
    if isnothing(product)
        product = (x,y) -> x.*y
    end
    ngrid = length(fongrid)
    #pairwisenorms = Matrix(undef,ngrid, ngrid)
    pairwisenorms = Matrix{typeof(product(fongrid[1], fongrid[1]))}(undef, (ngrid, ngrid))
    for i1 = 1:ngrid, i2 = 1:ngrid
        #for 
            #pairwisenorms[i1, i2] = transpose(fongrid[i1]) * weight * fongrid[i2]
            pairwisenorms[i1, i2] = product(fongrid[i1], fongrid[i2])
        #end
    end
    # pairwisenorms .= pairwisenorms .+ conj(permutedims(pairwisenorms, [2, 1]))# .- Diagonal(diag(pairwisenorms))
    # for ii = 1:ngrid
    #     pairwisenorms[ii, ii] = 0.5 .* pairwisenorms[ii, ii]
    # end
    return pairwisenorms
end

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
    n = length(terms[1][1])

    a = Vector{Int64}(undef,2)
    b = Vector{Int64}(undef,2)

    for (i, rowi) = enumerate(eachrow(terms)), (j, rowj) = enumerate(eachrow(terms))
        @inbounds if nonzero_rows[i] == false
            continue
        end
        #for (j, rowj) = enumerate(eachrow(terms))
            #@inbounds if iszero(pairwisenorms[maprowtouniquept[i], maprowtouniquept[j]])
            #    display(nonzero_rows[j])
            @inbounds if nonzero_rows[j] == false
                continue
            end
            #@inbounds @fastmath valstomultiply[i,j] = prod(productintegrals[rowi[1][k][1], rowj[1][k][1], rowi[1][k][2], rowj[1][k][2]] for k in eachindex(rowi[1]))
            # valstomultiply[1] = 1.0
            rowi1 = rowi[1]
            rowj1 = rowj[1]
            @inbounds for k in eachindex(rowi1)
                a .= rowi1[k]
                b .= rowj1[k]
                @fastmath valstomultiply[i,j] *= productintegrals[k][a[1],b[1],a[2],b[2]]
            end
            @inbounds feval .= feval .+ coeff_per_term[i] * coeff_per_term[j] * valstomultiply[i,j] * pairwisenorms[maprowtouniquept[i], maprowtouniquept[j]]
        #end
    end
    #@inbounds feval .= sum(sum(coeff_per_term[i] * coeff_per_term[j] * valstomultiply[i,j] * pairwisenorms[maprowtouniquept[i], maprowtouniquept[j]] for j=1:size(terms,1)) for i = 1:size(terms,1))

    feval = sqrt.(abs.(feval)) # abs prevents very small negative numbers (≈ 0) giving error

    if returnpairwisenorms == true
        return (feval, pairwisenorms)
    else
        return feval
    end
end

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

function computesparsegridc(I)
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

function generatebinaryvectors(n)
    vectors = Matrix{Integer}(undef, n, 0)
    for i in 0:2^n-1
        binary_rep = reverse(digits(i, base=2, pad=n))
        vectors = hcat(vectors, binary_rep)
    end
    return vectors
end

function createlagrangepolys(pts, space=Chebyshev(-1.0..1.0))
    p = Vector{Any}(undef, length(pts))
    for ii in eachindex(pts)
        # p[ii] = stablepoly(pts, ii)
        p[ii] = poly_Fun(pts, ii, space=space)
    end
    return p
end

fastvectorgenerator(k, n) = Iterators.product(repeat([1:(k+1)], n)...)

function createsmolyakmiset(n, k)
    iteratorI = fastvectorgenerator(k, n)
    accepted = Vector{NTuple{n,Int64}}()
    l = n + k
    for vector in iteratorI
        if sum(vector) <= l
            push!(accepted, vector)
        end
    end

    I = hcat(collect.(accepted)...)
    I = sortmi(I)
    return I
end

# function createsmolyakmiset(n,k)
#     vector = ones(n)
#     finished = false
#     accepted = Vector{Vector{Int}}()
#     l = n + k
#     mi_norm = n
#     while !finished
#         if mi_norm <= l
#             push!(accepted, vector)
#         end

#         vector[end] +=1
#         mi_norm = mi_norm + 1
#         if mi_norm >= l
#             for ii = n:-1:2
#                 if vector[ii] > l
#                     vector[ii] = 1
#                     vector[ii-1] += 1
#                     mi_norm = mi_norm - l + 1
#                 end
#             end
#         end
#         if vector[1] > l
#             finished = true
#         end
#     end
#     # return accepted = sortmi(hcat(collect.(accepted)...))
#     return accepted = hcat(accepted...)
# end

function createtensormiset(n, k)
    iteratorI = fastvectorgenerator(k, n)
    accepted = Vector{NTuple{n,Int64}}()
    l = n + k
    for vector in iteratorI
        if maximum(vector) <= l
            push!(accepted, vector)
        end
    end

    I = hcat(collect.(accepted)...)
    I = sortmi(I)
    return I
end

function is_leq(α::Vector, β::Vector)
    # Check if β ≤ α element-wise
    return all(β[i] ≤ α[i] for i in 1:length(α))
end

function downwards_closed_set(input_matrix::Matrix)
    num_indices, num_columns = size(input_matrix)
    output_matrix = copy(input_matrix)  # Initialize output matrix with input matrix

    # Iterate over each column (multi-index) in the input matrix
    for col_idx in 1:num_columns
        α = input_matrix[:, col_idx]

        # Generate all possible multi-indices β such that β ≤ α
        for β in Iterators.product((1:α[i] for i in 1:num_indices)...)
            β_vector = collect(β)
            if is_leq(α, β_vector) && !any(all(β_vector .== output_matrix[:, j]) for j in 1:size(output_matrix, 2))
                output_matrix = hcat(output_matrix, β_vector)
            end
        end
    end
    output_matrix = reduce(hcat, unique(eachcol(output_matrix)))

    return sortmi(output_matrix)
end

function spaceselector(domain)
    if length(domain) != 2
        Base.throw_checksize_error("Domain must be 2 values")
    end
    xmin = domain[1]
    xmax = domain[2]
    domain = xmin..xmax

    if isinf(domain.left) && isinf(domain.right)
        space = Hermite()
    elseif iszero(domain.left) && isinf(domain.right)
        space = Laguerre()
    else
        space = Chebyshev(domain)
    end
    return space
end

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

function lagrange_evaluation!(y, pts::Vector{<:Real}, ii, targetpoints)
    # Check uniqueness of points
    n = length(pts)
    pts = Float64.(pts)
    # Denominator
    denom = one(eltype(pts))
    @inbounds for j in 1:n
        if j != ii
            denom *= pts[ii] - pts[j]
        end
    end
    # Numerator
    numer = Ref{Float64}(0.0)
    difference = Ref{Float64}(0.0)
    @inbounds for k in eachindex(targetpoints)
        numer[] = 1.0
        x = targetpoints[k]
        @inbounds for j in 1:n
            if j != ii
                difference[] = pts[j]
                difference[] *= -1.0
                difference[] += x
                numer[] *= difference[]
            end
        end
        y[k] = numer[] / denom
    end
end