"""
    adaptive_sparsegrid(f, ndims; maxpts=100, proftol=1e-4, rule=doubling, knots=ccpoints)

Constructs an adaptive sparse grid for approximating the function `f` in `ndims` dimensions.

# Arguments
- `f`: Function to be approximated taking arguments x with length(x)=ndims.
- `ndims`: Dimension of domain.
- `maxpts`: (Optional) Maximum number of points to include in the sparse grid. Default is `100`.
- `proftol`: (Optional) Tolerance for profits. Default is `1e-4`.
- `rule`: (Optional) Level function(s) for sparse grid. Default is `doubling`.
- `knots`: (Optional) Knot function(s) for sparse grid. Default is `ccpoints`.

# Returns
- `sg`: Final sparse grid used to approximate `f`
- `f_on_z`: Evaluations of `f` on grid points in sparse grid `sg`
"""
function adaptive_sparsegrid(f, ndims; maxpts = 100, proftol=1e-4, rule = doubling, knots = ccpoints)
    MI = create_smolyak_miset(ndims,0)
    sg = create_sparsegrid(MI; rule=rule, knots=knots)

    Z = get_grid_points(sg)

    f_on_z = [f(z) for z in Z]
    kk=0

    maxmi = 5
    pcl = precompute_lagrange_integrals(maxmi, knots, rule)

    while true
        if get_n_grid_points(sg) > maxpts
            @info "Max grid points reached"
            break
        end
        kk+=1
        RM = get_reduced_margin(MI)
        MI_enhanced = add_mi(MI, RM)
        sg_enhanced = create_sparsegrid(MI_enhanced; rule=rule, knots=knots)

        # SOLVE
        f_on_z_enhanced = adaptive_solve(sg,sg_enhanced, f_on_z)

        # Update precompute_lagrange_integrals if necessary
        if maximum([maximum(α) for α in get_mi(get_mi_set(sg_enhanced))])  > maxmi
            maxmi += 3
            pcl = precompute_lagrange_integrals(maxmi, knots, rule)
        end

        # ESTIMATE
        # Interpolate to new grid
        p_α = adaptive_estimate(sg, sg_enhanced, f_on_z_enhanced, pcl)

        # MARK
        α_marked = adaptive_mark(RM, p_α, θ)

        # TERMINATE
        if adaptive_terminate(sg, p_α, maxpts, proftol)
            break
        end

        # REFINE
        sg, f_on_z = adaptive_refine(sg, α_marked, rule, knots)

        @info "Iteration: "*string(kk)*"    Number of points: "*string(get_n_grid_points(sg))*"    Max profit: "*string(p_max)
    end

    @info "Finished in "*string(kk)*" iterations"

    return (sg, f_on_z)
end

"""
    adaptive_solve(sg,sg_enhanced, f_on_z)

Computes function f on the new grid points in sg_enhanced

# Arguments
- `sg`: Sparse grid.
- `sg_enhanced`: Enhanced sparse grid.
- `f_on_z`: Function evaluations on the sparse grid `sg`.

# Returns
- Function evaluations on the enhanced sparse grid `sg_enhanced`.
"""
function adaptive_solve(sg,sg_enhanced, f_on_z)
    Z_enhanced = get_grid_points(sg_enhanced)
    sg_map = mapfromto(sg,sg_enhanced)
    f_on_z_enhanced = Vector{typeof(f_on_z[1])}(undef,get_n_grid_points(sg_enhanced))
    # Fill old evaluations
    f_on_z_enhanced[sg_map] = f_on_z
    # Fill new evaluations
    new_points = setdiff(1:get_n_grid_points(sg_enhanced), sg_map)
    f_on_z_enhanced[new_points] = [f(z) for z in Z_enhanced[new_points]]
    return f_on_z_enhanced
end

"""
    adaptive_estimate(sg, sg_enhanced, f_on_z_enhanced, pcl)

Estimates the profit of adding multi-indices {α} in reduced margin to the sparse grid.

# Arguments
- `sg`: Sparse grid.
- `sg_enhanced`: Enhanced sparse grid.
- `f_on_z_enhanced`: Function evaluations on the enhanced sparse grid.
- `pcl`: Precomputed Lagrange integrals.

# Returns
- Vector of profits for each multi-index α.
"""
function adaptive_estimate(sg, sg_enhanced, f_on_z_enhanced, pcl)
    Z_enhanced = get_grid_points(sg_enhanced)
    f_sg_on_z_enhanced = interpolate_on_sparsegrid(sg,f_on_z, Z_enhanced)
    p_α = Vector{Float64}(undef, length(get_mi(RM)))
    for (i,α) in enumerate(get_mi(RM))
        sg_α = create_sparsegrid(add_mi(MI,α), rule=rule, knots=knots)
        sg_map_α = mapfromto(sg_α,sg_enhanced)
        f_on_z_α = f_on_z_enhanced[sg_map_α]
        
        f_diff_α = f_on_z_α - f_sg_on_z_enhanced[sg_map_α]

        p_α[i] = compute_profit(sg_α, f_diff_α,pcl)
    end
    return p_α
end

"""
    compute_profit(sg_α, f_diff_α, pcl)

Computes the "profit" of a sparse grid supplemented by a multi-index α

# Arguments
- `sg_α`: Enhanced sparse grid with α
- `f_diff_α`: Function differences (surpluses) at the enhanced sparse grid points.
- `pcl`: Precomputed Lagrange integrals

# Returns
- Computed profit as Expected change in approximation.
"""
function compute_profit(sg_α, f_diff_α, pcl)
    return integrate_on_sparsegrid(sg_α, abs.(f_diff_α), pcl)
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
    α_marked = α[idx_sorted[1:k]]
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
function adaptive_refine(sg, α_marked, rule, knots)
    MI = get_mi_set(sg)
    MI = add_mi(MI, MISet(α_max))
    sg = create_sparsegrid(MI; rule=rule, knots=knots)
    sg_map_refine = mapfromto(sg,sg_enhanced)
    f_on_z = f_on_z_enhanced[sg_map_refine]
return sg, f_on_z