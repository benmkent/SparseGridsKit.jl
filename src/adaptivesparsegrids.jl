"""
    adaptive_sparsegrid(f, ndims; maxpts = 100, proftol = 1e-4, rule=Doubling(), knots=CCPoints(), θ=1e-4, type=:deltaint)

Constructs an adaptive sparse grid for approximating the function `f` in `ndims` dimensions.

# Arguments
- `f`: Function to be approximated taking arguments x with length(x)=ndims.
- `ndims`: Dimension of domain.
- `maxpts`: (Optional) Maximum number of points to include in the sparse grid. Default is `100`.
- `proftol`: (Optional) Tolerance for profits. Default is `1e-4`.
- `rule`: (Optional) Level function(s) for sparse grid. Default is `Doubling()`.
- `knots`: (Optional) Knot function(s) for sparse grid. Default is `CCPoints()`.
- `θ`: (Optional) Threshold for marking. Default is `1e-4`.
- `type`: (Optional) Type of profit computation. Default is `:deltaint`.
- `costfunction` (Optional) Cost function for fidelity related multi-indices.

# Returns
- `sg`: Final sparse grid used to approximate `f`
- `f_on_z`: Evaluations of `f` on grid points in sparse grid `sg`
"""
function adaptive_sparsegrid(f, ndims; maxpts = 100, proftol=1e-4, rule = Doubling(), knots = CCPoints(), θ=1e-4, type=:deltaint, costfunction=nothing)
    MI = create_smolyak_miset(ndims,0)
    sg = create_sparsegrid(MI; rule=rule, knots=knots)

    Z = get_grid_points(sg)

    f_on_z = [f(z) for z in Z]
    kk=0

    maxmi = 5
    pcl = precompute_lagrange_integrals(maxmi, knots, rule)

    # Set up data store for evaluation recycling
    datastore = EvaluationDictionary(sg, f_on_z)

    while true
        if get_n_grid_points(sg) > maxpts
            @info "Max grid points reached"
            break
        end
        kk+=1
        MI = get_mi_set(sg)
        RM = get_reduced_margin(MI)
        MI_enhanced = add_mi(MI, RM)
        sg_enhanced = create_sparsegrid(MI_enhanced; rule=rule, knots=knots)

        # SOLVE
        adaptive_solve!(f, datastore, sg_enhanced)

        # Update precompute_lagrange_integrals if necessary
        if maximum([maximum(α) for α in get_mi(get_mi_set(sg_enhanced))])  > maxmi
            maxmi += 3
            pcl = precompute_lagrange_integrals(maxmi, knots, rule)
        end

        # ESTIMATE
        p_α = adaptive_estimate(sg, datastore, pcl, rule, knots, type=type, costfunction=costfunction)

        # MARK
        α_marked = adaptive_mark(RM, p_α, θ)

        # TERMINATE
        if adaptive_terminate(sg, p_α, maxpts, proftol)
            break
        end

        # REFINE
        sg, f_on_z = adaptive_refine(sg, datastore, α_marked, rule, knots)

        @info "Iteration: "*string(kk)*"    Number of points: "*string(get_n_grid_points(sg))*"    Max profit: "*string(maximum(p_α))
    end

    @info "Finished in "*string(kk)*" iterations"

    return (sg, f_on_z)
end

"""
    adaptive_solve!(f, datastore, sg)

Adds function evaluations on sg to datastore

# Arguments
- `f`: Function
- `datastore`: EvaluationDictionary
- `sg`: Sparse grid

# Returns
- `datastore`: Updated EvaluationDictionary
"""
function adaptive_solve!(f, datastore, sg)
    # Retrieve evaluations
    f_on_z = retrieve_evaluations(datastore, sg)

    # Fill new evaluations
    for (ii,z) in enumerate(get_grid_points(sg))
        if isnothing(f_on_z[ii])
            f_on_z[ii] = f(z)
            add_evaluation!(datastore, z, f_on_z[ii])
        end
    end
end

"""
    adaptive_estimate(sg, datastore, pcl, type=:deltaint)

Estimates the profit of adding multi-indices {α} in reduced margin to the sparse grid.

# Arguments
- `sg`: Sparse grid.
- `datastore`: EvaluationDictionary.
- `pcl`: Precomputed Lagrange integrals.
- `rule`: Level function(s) for sparse grid.
- `knots`: Knot function(s) for sparse grid.
- `type`: Type of profit computation. Default is `:deltaint`.

# Returns
- Vector of profits for each multi-index α.
"""
function adaptive_estimate(sg, datastore, pcl, rule, knots; type=:deltaint, costfunction=nothing)
    f_on_z = retrieve_evaluations(datastore, sg)

    MI = get_mi_set(sg)
    RM = get_reduced_margin(MI)
    
    p_α = Vector{Float64}(undef, length(get_mi(RM)))
    for (i,α) in enumerate(get_mi(RM))
        sg_α = create_sparsegrid(add_mi(MI,α), rule=rule, knots=knots)
        f_on_z_α = retrieve_evaluations(datastore, sg_α)
        
        # cost = get_n_grid_points(sg_α) - get_n_grid_points(sg)
        cost = length(setdiff(get_grid_points(sg_α), get_grid_points(sg)))
        if costfunction != nothing
            α_fidelities = α[isa.(knots,FidelityPoints)]
            costpersolve = costfunction(α_fidelities)
            cost = cost*costpersolve
        end

        p_α[i] = compute_profit(sg_α, sg, f_on_z_α, f_on_z, cost, pcl, type=type)
    end
    return p_α
end

"""
    compute_profit(sg_α, f_α, f, cost, pcl; type=:deltaint)

Computes the "profit" of a sparse grid supplemented by a multi-index α

# Arguments
- `sg_α`: Enhanced sparse grid with α
- `sg`: Sparse grid
- `f_α`: Function evaluations at grid points in α
- `f`: Function evaluations on the sparse grid.
- `cost`: Cost of adding α to the sparse grid.
- `pcl`: Precomputed Lagrange integrals
- `type`: Type of profit computation. Default is `:deltaint`. Options are `:deltaint`, `:deltaintcost`, `:Linf`, `:Linfcost`.

# Returns
- Computed profit as Expected change in approximation.
"""
function compute_profit(sg_α, sg, f_α, f, cost, pcl; type=:deltaint)
    if type == :deltaint
        profit =  abs(integrate_on_sparsegrid(sg_α, f_α, pcl) - integrate_on_sparsegrid(sg, f, pcl))
    elseif type == :deltaintcost
        profit =  abs(integrate_on_sparsegrid(sg_α, f_α, pcl) - integrate_on_sparsegrid(sg, f, pcl))/cost
    elseif type == :Linf
        f_interp = interpolate_on_sparsegrid(sg, f, get_grid_points(sg_α))
        profit = maximum(abs.(f_α - f_interp))
    elseif type == :Linfcost
        f_interp = interpolate_on_sparsegrid(sg, f, get_grid_points(sg_α))
        profit = maximum(abs.(f_α - f_interp))/cost
    else
        @warn "Invalid profit type " * string(type) * ". Using :deltaint"
        profit =  abs(integrate_on_sparsegrid(sg_α, f_α, pcl) - integrate_on_sparsegrid(sg, f, pcl))
    end
    return profit
end

"""
    adaptive_mark(α, p_α, θ=1e-4)

Dorfler style marking of multi-indices based on profits

# Arguments
- `α`: Multi-indices in reduced margin.
- `p_α`: Vector of profits.
- `θ`: Threshold for marking. Default is `1e-4`.
"""
function adaptive_mark(α, p_α, θ=1e-4)
    idx_sorted = sortperm(p_α, rev=true) 
    # Sort indices by decreasing p_α
    p_sorted = p_α[idx_sorted]
        
    p_total = sum(p_sorted)
    p_cumsum = cumsum(p_sorted)
    
    # Dorfler marking
    k = findfirst(p_cumsum .≥ θ * p_total)
    α_marked = MISet(get_mi(α)[idx_sorted[1:k]])
    return α_marked
end

"""
   adaptive_terminate(sg, p_α, maxpts, proftol)

Test profits and sparse grid to determine loop termination

# Arguments
- `sg`: Sparse grid.
- `p_α`: Vector of profits.
- `maxpts`: Maximum number of grid points allowed.
- `proftol`: Profit tolerance.

# Returns
- Flag to indicate loop termination if `maxpts` has been reached, or `proftol` attained.
"""
function adaptive_terminate(sg, p_α, maxpts, proftol)
    retval = false
    if get_n_grid_points(sg) >= maxpts
        @info "Reached max points: "*string(maxpts)
        retval = true
    elseif all(proftol .> p_α)
        @info "Profits below proftol: "*string(proftol)
        retval= true
    end
    return retval
end

"""
    adaptive_refine(sg, α_marked, rule, knots)

Refines the sparse grid based on marked multi-indices

# Arguments
- `sg`: Sparse grid.
- `α_marked`: Marked multi-indices.
- `rule`: Level function(s) for sparse grid.
- `knots`: Knot function(s) for sparse grid.

# Returns
- `sg`: Refined sparse grid.
- `f_on_z`: Function evaluations on the refined sparse grid.
"""
function adaptive_refine(sg, datastore, α_marked, rule, knots)
    MI = get_mi_set(sg)
    MI = add_mi(MI, α_marked)
    sg = create_sparsegrid(MI; rule=rule, knots=knots)
    f_on_z = retrieve_evaluations(datastore, sg)
    return sg, f_on_z
end

"""
    evaluations_database
"""
mutable struct EvaluationDictionary
    dictionary
end


"""
    EvaluationDictionary()
    
Creates an empty EvaluationDictionary

# Returns
- `evaluations`: EvaluationDictionary
"""
function EvaluationDictionary()
    return EvaluationDictionary(Dict())
end

"""
    EvaluationDictionary(sg,f_on_z)

Creates an EvaluationDictionary from a sparse grid and function evaluations

# Arguments
- `sg`: Sparse grid
- `f_on_z`: Function evaluations on the sparse grid

# Returns
- `evaluations`: EvaluationDictionary
"""
function EvaluationDictionary(sg,f_on_z)
    evaluations = EvaluationDictionary()
    add_evaluations!(evaluations, sg, f_on_z)
    return evaluations
end

"""
    add_evaluations!(evaluations, sg::SparseGrid, f_on_z)

Adds the evaluations f_on_z into the evaluations_database

# Arguments
- `evaluations`: EvaluationDictionary
- `sg`: Sparse grid
- `f_on_z`: Function evaluations on the sparse grid

# Returns
- `evaluations`: Updated EvaluationDictionary
"""
function add_evaluations!(evaluations, sg::SparseGrid, f_on_z)
    for (i, z) in enumerate(get_grid_points(sg))
        if !haskey(evaluations.dictionary, z)
            evaluations.dictionary[z] = f_on_z[i]
        end
    end
    return evaluations
end

"""
    add_evaluation!(evaluations, z, f_at_z)

Adds the evaluation f_at_z into the evaluations_database

# Arguments
- `evaluations`: EvaluationDictionary
- `z`: Point
- `f_at_z`: Function evaluation at z

# Returns
- `evaluations`: Updated EvaluationDictionary
"""
function add_evaluation!(evaluations, z, f_at_z)
    if !haskey(evaluations.dictionary, z)
        evaluations.dictionary[z] = f_at_z
    end
    return evaluations
end

"""
    retrieve_evaluations(evaluations, sg)

Retrieves the evaluations from the evaluations_database

# Arguments
- `evaluations`: EvaluationDictionary
- `sg`: Sparse grid

# Returns
- `f_on_z`: Function evaluations on the sparse grid
"""
function retrieve_evaluations(evaluations, sg)
    f_on_z = Vector(undef, get_n_grid_points(sg))
    for (i, z) in enumerate(get_grid_points(sg))
        if haskey(evaluations.dictionary, z)
            f_on_z[i] = evaluations.dictionary[z]
        else
            f_on_z[i] = nothing
        end
    end
    return f_on_z
end