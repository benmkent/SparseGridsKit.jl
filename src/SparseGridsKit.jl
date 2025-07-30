module SparseGridsKit

export SparseGrid
export MISet
export SparseGridApproximation
export SpectralSparseGridApproximation

export create_sparsegrid, get_grid_points, get_n_grid_points, map_from_to, get_mi_set
export add_mi, get_mi, get_n_mi, get_margin, get_reduced_margin, create_tensor_miset, create_smolyak_miset, create_box_miset, create_totaldegree_miset, create_rule_miset
export check_admissibility, check_index_admissibility
export interpolate_on_sparsegrid, integrate_on_sparsegrid, integrate_L2_on_sparsegrid, compute_quadrature_weights!
export precompute_lagrange_integrals, precompute_pairwise_norms
export adaptive_sparsegrid
export convert_to_spectral_approximation
export derivative

#export ccpoints, uniformpoints, lejapoints, transformdomain, gausshermitepoints, gausslegendrepoints, lejapoints
export CCPoints, UniformPoints, GaussLegendrePoints, GaussHermitePoints, LejaPoints, CustomPoints, Points
export Doubling, Linear, Tripling, TwoStep, CustomLevel, Level

export Fidelity, FidelityPoints, multifidelityfunctionwrapper

export genz

## Use functionality defined in import file
include("sparsegrids.jl")
include("genz.jl")
include("misets.jl")
include("knots.jl")
include("knots_structures.jl")
include("adaptivesparsegrids.jl")
include("spectralsparsegrids.jl")
include("verifyinputs.jl")
include("multifidelity.jl")
include("differentiation.jl")
# Include plots last
include("sparsegridplots.jl")

"""
    precompute_lagrange_integrals(max_mi)

Precomputes the product integrals for Lagrange basis functions up to a given maximum multi-index (`max_mi`).

# Arguments
- `max_mi`: The maximum multi-index for which to precompute integrals.

# Returns
- A vector of precomputed product integrals for the Lagrange basis.
"""
function precompute_lagrange_integrals(max_mi, knots=CCPoints(), rule=Doubling())
    precompute = sparsegridprecompute(max_mi, knots, rule)
    return precompute.productintegrals
end

"""
    create_sparsegrid(mi_set; rule=Doubling(), knots=CCPoints())

Creates a sparse grid based on the provided multi-index set (`mi_set`).

# Arguments
- `mi_set`: An instance of `MISet` containing the multi-index set for grid construction.
- `rule`: Map from index to number of points. Optional (default Doubling()). Function or vector of functions (for each dimension).
- `knots`: Type of knots. Optional (default Clenshaw--Curtis). Function or vector of functions (for each dimension).

# Returns
- A sparse grid object constructed using the specified multi-index set.
"""
function create_sparsegrid(mi_set; rule=Doubling(), knots=CCPoints())
    verifyinputs!(mi_set)
    #verifyinputs(knots)

    miset_matrix = hcat(mi_set.mi...)
    sg = createsparsegrid(miset_matrix; rule=rule, knots=knots)
    return sg
end

"""
    map_from_to(sg_from, sg_to)

Maps data from one sparse grid (`sg_from`) to another (`sg_to`).

# Arguments
- `sg_from`: The source sparse grid.
- `sg_to`: The target sparse grid.

# Returns
- A vector that maps data from `sg_from` to `sg_to`.
"""
function map_from_to(sg_from, sg_to)
    map_vector = mapfromto(sg_from, sg_to)
    return map_vector
end

"""
    get_grid_points(sg)

Retrieves the grid points from a sparse grid (`sg`).

# Arguments
- `sg`: A sparse grid object.

# Returns
- A vector of vectors, where each inner vector represents a grid point.
"""
function get_grid_points(sg)
    matrix_points = sg.grid_points
    grid_point_vector = [Vector(v) for v in eachrow(matrix_points)]
    @debug "Returned "*string(length(grid_point_vector))*" sparse grid points"
    return grid_point_vector
end

"""
    get_n_grid_points(sg)

Returns the number of grid points in the sparse grid (`sg`).

# Arguments
- `sg`: A sparse grid object.

# Returns
- The number of grid points in the sparse grid.
"""
function get_n_grid_points(sg)
    n = size(sg.grid_points, 1)
    @debug "Number of grid points "*string(n)
    return n
end

"""
    get_mi_set(sg)

Generates a downwards-closed set of multi-indices from a sparse grid (`sg`).

# Arguments
- `sg`: A sparse grid object.

# Returns
- An `MISet` containing the downwards-closed set of multi-indices.
"""
function get_mi_set(sg)
    mi = sg.multi_index_set
    mi_matrix = mi
    @debug "Getting MI set, currently has "*string(size(mi_matrix),1)*" terms"
    mi_matrix_dc = downwards_closed_set(mi_matrix)
    @debug "Enforced downwards close, now  "*string(size(mi_matrix_dc,1))*" terms"
    mi_set = MISet([Vector(v) for v in eachcol(mi_matrix_dc)])
    return mi_set
end

"""
    interpolate_on_sparsegrid(sg, f_on_grid, target_points)

Interpolates a function (`f_on_grid`) defined on a sparse grid (`sg`) to a set of target points.

# Arguments
- `sg`: The source sparse grid.
- `f_on_grid`: A vector of function values on the sparse grid.
- `target_points`: A vector of target points for interpolation.

# Returns
- A vector of interpolated values at the target points.
"""
function interpolate_on_sparsegrid(sg, f_on_grid, target_points)
    #target_points_matrix = hcat(target_points...)'

    flagmatrix=false

    if !(target_points isa AbstractVector{<:AbstractVector{<:Real}})
        if target_points isa AbstractVector{<:Real}
            target_points = [target_points]
        elseif target_points isa Real
            target_points = [[target_points]]
        elseif target_points isa AbstractMatrix{<:AbstractVector{<:Real}}
            flagmatrix = true
            matrix_shape = size(target_points)
            target_points = vec(target_points)
        else
            error("target_points must be a vector of vectors, a vector of reals, or a matrix of vectors.")
        end
    end

    data_converted = false
    if f_on_grid isa AbstractVector{<:Real} && length(f_on_grid) == get_n_grid_points(sg)
        # Data converted to vector of vectors
        data_converted = true
        f_on_grid = [[fx] for fx in f_on_grid]
    end

    verifyinputs(sg.domain, target_points)

    @debug "Interpolating onto "*string(length(target_points))*" target points"
    f_on_target_points = interpolateonsparsegrid(sg, f_on_grid, target_points)

    if flagmatrix
        f_on_target_points = reshape(f_on_target_points, matrix_shape)
    end
    if data_converted
        f_on_target_points = [fx[1] for fx in f_on_target_points]
    end
    return f_on_target_points
end
"""
    interpolate_on_sparsegrid(sga::SparseGridApproximation, target_points)

Interpolates a SparseGridApproximation to a set of target points.

# Arguments
- `sga`: SparseGridApproximation.
- `target_points`: A vector of target points for interpolation.

# Returns
- A vector of interpolated values at the target points.
"""
function interpolate_on_sparsegrid(sga::SparseGridApproximation, target_points)
    interpolate_on_sparsegrid(sga.sg, sga.fongrid, target_points)
end

"""
    integrate_on_sparsegrid(sg, f_on_grid, precomputed_lagrange_integrals)

Integrates a function (`f_on_grid`) over a sparse grid (`sg`) using precomputed Lagrange integrals.

# Arguments
- `sg`: The sparse grid for integration.
- `f_on_grid`: A vector of function values on the sparse grid.

# Returns
- The integral of the function over the sparse grid.
"""
function integrate_on_sparsegrid(sg, f_on_grid)
    if isempty(sg.quadrature_weights)
        compute_quadrature_weights!(sg)
    end
    integral = sum(sg.quadrature_weights .* f_on_grid)
    return integral
end

""" compute_quadrature_weights(sg)
Computes the quadrature weights for a sparse grid (`sg`).
# Arguments
- `sg`: The sparse grid for which to compute the quadrature weights.
# Returns
- Updates the `quadrature_weights` field of the sparse grid with computed weights.
"""
# function compute_quadrature_weights!(sg)
#     # Compute quadrature weight for each term
#     q = ones(length(sg.data.terms))
#     weights = zeros(get_n_grid_points(sg))
#     for (ii,t) in enumerate(sg.data.terms)
#         for (jj,tii) in enumerate(t)
#             q[ii] *= sg.data.weight_sequences_per_dimension[jj][tii[1]][tii[2]]
#         end
#         weights[sg.data.terms_to_grid_points[ii]] += q[ii] * sg.data.coeff_per_term[ii]
#     end
#     sg.quadrature_weights = weights
# end

function compute_quadrature_weights!(sg)
    weights = zeros(get_n_grid_points(sg))

    @inbounds for ii in 1:length(sg.data.terms)
        t = sg.data.terms[ii]
        qi = 1.0
        for jj in 1:length(t)
            wi = sg.data.weight_sequences_per_dimension[jj]
            qi *= wi[t[jj][1]][t[jj][2]]
        end
        weights[sg.data.terms_to_grid_points[ii]] += qi * sg.data.coeff_per_term[ii]
    end

    sg.quadrature_weights = weights
end


"""
    integrate_L2_on_sparsegrid(sg, f_on_grid, precomputed_lagrange_integrals; product=dot, precomputed_pairwise_norms=nothing)

Computes the L2 norm of a function (`f_on_grid`) over a sparse grid (`sg`) using precomputed integrals.

# Arguments
- `sg`: The sparse grid.
- `f_on_grid`: A vector of function values on the sparse grid.
- `precomputed_lagrange_integrals`: Precomputed Lagrange integrals.
- `product`: A function to compute the inner product (default is `dot`).
- `precomputed_pairwise_norms`: Optional precomputed norms for optimization.

# Returns
- The L2 norm of the function over the sparse grid.
"""
function integrate_L2_on_sparsegrid(sg, f_on_grid, precomputed_lagrange_integrals; product=nothing, precomputed_pairwise_norms=nothing)
    precompute = (;productintegrals=precomputed_lagrange_integrals)
    integral = L2onsparsegrid(sg, f_on_grid, precompute; product=product, pairwisenorms=precomputed_pairwise_norms)
    return integral
end

"""
    precompute_pairwise_norms(f_on_grid, product)

Precomputes pairwise norms for function values (`f_on_grid`) using a specified product function.

# Arguments
- `f_on_grid`: A vector of function values on the sparse grid.
- `product`: A function to compute the product between pairs of function values.

# Returns
- A matrix of pairwise norms.
"""
function precompute_pairwise_norms(f_on_grid; product=nothing)
    pairwisenorms = computepairwisenorms(f_on_grid, product=product)
    return pairwisenorms
end

end