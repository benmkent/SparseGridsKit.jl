using LinearAlgebra
using StaticArrays
using Polynomials

import Base.:+
import Base.:-

# mutable struct SparseTerm
#     data::Vector{SVector{2, Int}}
# end

struct SparseTerm{N}
    data::NTuple{N, SVector{2,Int}}
end

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

function SparseGrid(dims, domain, multi_index_set, combination_coeff, grid_points, knots, rules, data)
    sg = SparseGrid(dims, domain, multi_index_set, combination_coeff, grid_points, [], knots, rules, data)
    compute_quadrature_weights!(sg)
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

    data = sparse_grid_data{dims}(terms, points_to_unique_indices, terms_to_grid_points, coeff_per_term, terms_to_unique_points_per_dimension, sg.data.point_sequences_per_dimension, sg.data.unique_points_per_dimension, weight_sequences_per_dimension)

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
        (unique_points_per_dimension, point_sequences_per_dimension, terms_to_unique_points_per_dimension, weight_sequences_per_dimension) = sparsegridonedimdata(maxmi, knots[ii], rule[ii], knots[ii].domain)
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
    (fromptids, _, _) = sparsegridpoints(sgfrom.data.terms, sgto.data.terms_to_unique_points_per_dimension, sgto.data.unique_points_per_dimension, sgto.dims)

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

# function sparsegridpoints(grid, terms_to_unique_points_per_dimension::Vector, unique_points_per_dimension, dims)
#     points_data = sparsegridpoints(grid, terms_to_unique_points_per_dimension, unique_points_per_dimension, dims, false)
#     return points_data
# end

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
    # points_to_unique_indices = reduce(vcat, [collect(t)' for t in points_to_unique_indices])
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

function sparsegridterms(multi_index_set, point_sequences_per_dimension; combination_coeff=nothing)
    if isnothing(combination_coeff)
        combination_coeff = computesparsegridc(multi_index_set)
    end
    n = size(multi_index_set,1)
    # terms = Vector{SVector{n,Int}}(undef, 0)
    # terms = Vector{NTuple{n,Vector{Int}}}(undef, 0)
    nterms_batch = 1000
    terms = Vector{SparseTerm{n}}(undef,nterms_batch)
    coeff_per_term = Vector{Int}(undef,nterms_batch)

    count = 1
    length_terms = nterms_batch

    for (miind, mi) in enumerate(eachcol(multi_index_set))
        if true #combination_coeff[miind] != 0
            # indexarray = Vector{Vector{Vector{Int}}}(undef, 0)
            # for (ind,i) in enumerate(mi)
            #     indexarrayi = Vector{Vector{Int}}(undef, 0)
            #     if isa(point_sequences_per_dimension[1][1],Number)
            #         for j = eachindex(point_sequences_per_dimension[i])
            #             push!(indexarrayi, [i, j])
            #         end
            #     else
            #         for j = eachindex(point_sequences_per_dimension[ind][i])
            #             push!(indexarrayi, [i, j])
            #         end
            #     end
            #     push!(indexarray, collect(indexarrayi))
            # end
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
                # push!(terms, NTuple{n, Vector{Int}}(Tuple(row)))
                # row_svec = {SVector{2,Int}}(undef, n)
                # for (i,r) in enumerate(row)
                #     row_svec[i] = SVector{2,Int}(r[1],r[2])
                # end
                # push!(terms, SparseTerm{n}(row_svec))
                # push!(terms, SparseTerm([SVector{2,Int}(r[1],r[2]) for r in row]))
                # push!(terms, SparseTerm{n}(ntuple(i -> SVector{2,Int}(row[i][1], row[i][2]), n)))
                # nt = ntuple(i -> SVector{2,Int}(row[i][1], row[i][2]), n)
                # t = SparseTerm{n}(nt)
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
    # Clean up multi_index_set and combination_coeff
    #zeroindices = combination_coeff .!= 0
    #multi_index_set = multi_index_set[:, zeroindices]
    #combination_coeff = combination_coeff[zeroindices]
   

    return (terms, coeff_per_term, combination_coeff, multi_index_set)
end

@generated function build_sparse_term(::Val{n}, row) where n
    # Unroll tuple construction for performance
    return :(SparseTerm{$n}(tuple($(map(i -> :(SVector{2,Int}(row[$i][1], row[$i][2])), 1:n)...))))
end

@noinline function sparsegridonedimdata(maxmi, knots, rule, domain)
    unique_points_per_dimension = Vector{Float64}(undef, 0)
    point_sequences_per_dimension = Vector{Vector{Float64}}(undef, maxmi)
    weight_sequences_per_dimension = Vector{Vector{Float64}}(undef, maxmi)
    #poly = Vector{Vector{Polynomial{Float64,:x}}}(undef, maxmi)
    terms_to_unique_points_per_dimension = Dict{Vector{Int},Int}()

    # Based on knots define appropriate ApproxFun space for underlying polynomials
    # space = spaceselector(domain)

    for ii = 1:maxmi
        # Detect special case of fidelity
        if isa(knots, FidelityPoints)
            point_sequences_per_dimension[ii], weight_sequences_per_dimension[ii] = knots(ii)
        else
            point_sequences_per_dimension[ii], weight_sequences_per_dimension[ii] = knots(rule(ii))
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
    return (unique_points_per_dimension, point_sequences_per_dimension, terms_to_unique_points_per_dimension, weight_sequences_per_dimension)
end

function interpolateonsparsegrid(sparsegrid, fongrid, targetpoints; evaltype=nothing)  
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

    # Precompute lagrange interpolation once
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

    #for k = 1:length(feval)
        for i in eachindex(terms)
        #targetpt = targetpoints[k]
        #for (i, row) = enumerate(eachrow(terms))
            #val = c[i]*prod([polyperlevel[row[1][j][1]][row[1][j][2]](newpt[1][j]) for j in eachindex(row[1])])
            # @inbounds  cweightedpolyprod[k,i] = coeff_per_term[i] * prod(poly[j][row[1][j][1]][row[1][j][2]](targetpt[j]) for j in eachindex(row[1]))
            # row = @view terms[i,:]
            row = terms[i].data
            @inbounds begin
                acc .= 1.0
                for j in eachindex(row)
                    # f = poly[j][row[1][j][1]][row[1][j][2]]
                    # f_placeholder = f.([targetpoints[k][j] for k = 1:ntarget])
                    # targetpoints_j = [targetpoints[k][j] for k = 1:ntarget]
                    @inbounds targetpoints_j = @view targetpoints_matrix[:,j]
                    r1 = row[j][1]
                    r2 = row[j][2]
                    pts = sparsegrid.data.point_sequences_per_dimension[j][r1]
                    
                    # lagrange_evaluation!(f_placeholder,pts,r2, targetpoints_j)
                    # display(norm(f_placeholder - lagrange_evaluation_precompute[j][r1][r2]))
                    f_placeholder = lagrange_evaluation_precompute[j][r1][r2]
                    acc .*= f_placeholder
                end
                # cweightedpolyprod[:,i] = coeff_per_term[i] * acc
                mul!(cweightedpolyprod_i,coeff_per_term[i],acc)
            end
            # feval += cweightedpolyprod[:,i] * fongrid[maprowtouniquept[i]]
            # if using_floats
                # for kk = eachindex(feval)
                #     temp_mul[] = fongrid[maprowtouniquept[i]]
                #     temp_mul[] = temp_mul[] * cweightedpolyprod_i[kk]
                #     @inbounds feval[kk] = feval[kk] + temp_mul[]
                # end
            # for kk = eachindex(feval)
            #     @inbounds copyto!(temp_mul[], fongrid[maprowtouniquept[i]])
            #     @inbounds for j in eachindex(temp_mul[])
            #         temp_mul[][j] *= cweightedpolyprod_i[kk]
            #     end
            #     @inbounds feval[kk] .+= temp_mul[]
            # end
            fused_update!(feval, fongrid[maprowtouniquept[i]], cweightedpolyprod_i,process_elementwise)
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
    n = length(terms[1].data)

    a = Vector{Int64}(undef,2)
    b = Vector{Int64}(undef,2)

    for (i, rowi) = enumerate(terms), (j, rowj) = enumerate(terms)
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
            rowi1 = rowi.data
            rowj1 = rowj.data
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
    (adm, mi_missing) = check_admissibility(MISet([vec(collect(c)) for c in eachcol(I)]))
    if adm
        coeff = computesparsegridc_admissible(I)
    else
        # Slower
        coeff = computesparsegridc_old(I)
    end
end

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

function generatebinaryvectors(n)
    vectors = Matrix{Integer}(undef, n, 0)
    for i in 0:2^n-1
        binary_rep = reverse(digits(i, base=2, pad=n))
        vectors = hcat(vectors, binary_rep)
    end
    return vectors
end

fastvectorgenerator(k, n) = Iterators.product(repeat([1:(k+1)], n)...)

function createsmolyakmiset(dim::Int, level::Int)
    multi_index_list = Vector{Vector{Int}}()
    current_index = Vector{Int}(undef, dim)

    total_order_min = dim
    total_order_max = dim + level

    for total_order in total_order_min:total_order_max
        generate_compositions!(multi_index_list, current_index, total_order, dim, 1)
    end

    return sortmi(hcat(multi_index_list...))
end

function generate_compositions!(
    output_list::Vector{Vector{Int}},
    composition_buffer::Vector{Int},
    remaining_sum::Int,
    total_parts::Int,
    current_position::Int,
)
    if current_position == total_parts
        composition_buffer[current_position] = remaining_sum
        push!(output_list, copy(composition_buffer))
        return
    end

    max_value = remaining_sum - (total_parts - current_position)
    for value in 1:max_value
        composition_buffer[current_position] = value
        generate_compositions!(
            output_list,
            composition_buffer,
            remaining_sum - value,
            total_parts,
            current_position + 1,
        )
    end
end

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

function is_leq(α::Vector, β::Vector)
    # Check if β ≤ α element-wise
    return all(β[i] ≤ α[i] for i in 1:length(α))
end

function downwards_closed_set(input_matrix::Matrix{Int})
    num_indices, num_columns = size(input_matrix)

    # Use a Set to store seen multi-indices as NTuples (cheaper lookup)
    seen = Set{NTuple{typeof(input_matrix[1,1]), Int}}()
    result = NTuple{num_indices, Int}[]

    for col_idx in 1:num_columns
        α = view(input_matrix, :, col_idx)

        # Generate all β ≤ α (componentwise)
        for β in Iterators.product((0:α[i] for i in 1:num_indices)...)
            β_tuple = β  # already an NTuple

            # Only keep β if β ≤ α (can omit if product bounds are strict)
            # if is_leq(α, β_tuple)   # Optional: may not be needed
            if !in(β_tuple, seen)
                push!(result, β_tuple)
                push!(seen, β_tuple)
            end
            # end
        end
    end

    # Convert result to a matrix
    output_matrix = hcat(result...)  # this is d × n
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
