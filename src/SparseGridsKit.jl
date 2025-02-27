module SparseGridsKit

export SparseGrid
export MISet
export SparseGridApproximation
export SpectralSparseGridApproximation

export create_sparsegrid, get_grid_points, get_n_grid_points, map_from_to, get_mi_set
export add_mi, get_mi, get_n_mi, get_margin, get_reduced_margin, create_tensor_miset, create_smolyak_miset, create_box_miset, create_totaldegree_miset, create_rule_miset
export interpolate_on_sparsegrid, integrate_on_sparsegrid, integrate_L2_on_sparsegrid
export precompute_lagrange_integrals, precompute_pairwise_norms
export adaptive_sparsegrid
export convert_to_spectral_approximation

export ccpoints, uniformpoints, linearpoints, lejapoints, transformdomain, gausshermitepoints, gausslegendrepoints
export doubling, linear, tripling, twostep

export genz

## Use functionality defined in import file
include("sparsegrids.jl")
include("genz.jl")
include("sparsegridplots.jl")
include("misets.jl")
include("knots.jl")
include("adaptivesparsegrids.jl")
include("spectralsparsegrids.jl")

"""
    precompute_lagrange_integrals(max_mi, domain)

Precomputes the product integrals for Lagrange basis functions up to a given maximum multi-index (`max_mi`).

# Arguments
- `max_mi`: The maximum multi-index for which to precompute integrals.

# Returns
- A vector of precomputed product integrals for the Lagrange basis.
"""
function precompute_lagrange_integrals(max_mi, domain, knots=ccpoints, rule=doubling)
    precompute = sparsegridprecompute(max_mi, domain, knots, rule)
    return precompute.productintegrals
end

"""
    create_sparsegrid(mi_set, domain)

Creates a sparse grid based on the provided multi-index set (`mi_set`).

# Arguments
- `mi_set`: An instance of `MISet` containing the multi-index set for grid construction.
- `rule`: Map from index to number of points. Optional (default doubling). Function or vector of functions (for each dimension).
- `knots`: Type of knots. Optional (default Clenshaw--Curtis). Function or vector of functions (for each dimension).

# Returns
- A sparse grid object constructed using the specified multi-index set.
"""
function create_sparsegrid(mi_set, domain; rule=doubling, knots=ccpoints)
    miset_matrix = hcat(mi_set.mi...)
    sg = createsparsegrid(miset_matrix, domain; rule=rule, knots=knots)
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
    matrix_points = sg.sparsegridpts
    grid_point_vector = [Vector(v) for v in eachrow(matrix_points)]
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
    return size(sg.sparsegridpts, 1)
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
    mi = sg.MI
    mi_matrix = mi
    mi_matrix_dc = downwards_closed_set(mi_matrix)
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
    f_on_target_points = interpolateonsparsegrid(sg, f_on_grid, target_points)
    return f_on_target_points
end

"""
    integrate_on_sparsegrid(sg, f_on_grid, precomputed_lagrange_integrals)

Integrates a function (`f_on_grid`) over a sparse grid (`sg`) using precomputed Lagrange integrals.

# Arguments
- `sg`: The sparse grid for integration.
- `f_on_grid`: A vector of function values on the sparse grid.
- `precomputed_lagrange_integrals`: A vector of precomputed Lagrange integrals.

# Returns
- The integral of the function over the sparse grid.
"""
function integrate_on_sparsegrid(sg, f_on_grid, precomputed_lagrange_integrals)
    precompute = (;productintegrals=precomputed_lagrange_integrals)
    integral = integrateonsparsegrid(sg, f_on_grid, precompute; evaltype=nothing)
    return integral
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