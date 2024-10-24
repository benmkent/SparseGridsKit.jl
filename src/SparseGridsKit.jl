module SparseGridsKit

export SparseGrid
export MISet

export create_sparsegrid, get_grid_points, get_n_grid_points, map_from_to, get_mi_set
export add_mi, get_mi, get_n_mi, get_margin, get_reduced_margin, create_smolyak_miset
export interpolate_on_sparsegrid, integrate_on_sparsegrid, integrate_L2_on_sparsegrid
export precompute_lagrange_integrals, precompute_pairwise_norms

## Use functionality defined in import file
include("sparsegrids.jl")

"""
    precompute_lagrange_integrals(max_mi)

Precomputes the product integrals for Lagrange basis functions up to a given maximum multi-index (`max_mi`).

# Arguments
- `max_mi`: The maximum multi-index for which to precompute integrals.

# Returns
- A vector of precomputed product integrals for the Lagrange basis.
"""
function precompute_lagrange_integrals(max_mi)
    precompute = sparsegridprecompute(max_mi)
    return precompute.productintegrals
end

"""
    create_sparsegrid(mi_set)

Creates a sparse grid based on the provided multi-index set (`mi_set`).

# Arguments
- `mi_set`: An instance of `MISet` containing the multi-index set for grid construction.

# Returns
- A sparse grid object constructed using the specified multi-index set.
"""
function create_sparsegrid(mi_set)
    miset_matrix = hcat(mi_set.mi...)
    sg = createsparsegrid(miset_matrix; rule=doubling)
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

mutable struct MISet
    mi::Vector
end

"""
    add_mi(mi_set, mi_set_new)

Adds a new set of multi-indices (`mi_set_new`) to an existing multi-index set (`mi_set`).

# Arguments
- `mi_set`: The original multi-index set.
- `mi_set_new`: A new `MISet` to be added to the original.

# Returns
- A new `MISet` containing the combined and sorted multi-indices.
"""
function add_mi(mi_set::MISet, mi_set_new::MISet)
    mi = copy(get_mi(mi_set))
    mi_new = get_mi(mi_set_new)
    append!(mi, mi_new)
    col_matrix = hcat(mi...)
    col_matrix_unique = unique(eachcol(col_matrix))
    col_matrix_sorted = sortmi(hcat(col_matrix_unique...))
    mi_set_new = MISet([Vector(v) for v in eachcol(col_matrix_sorted)])
    return mi_set_new
end

"""
    add_mi(mi_set, mi_new)

Adds a single new multi-index (`mi_new`) to an existing multi-index set (`mi_set`).

# Arguments
- `mi_set`: The original multi-index set.
- `mi_new`: A new multi-index vector to be added.

# Returns
- A new `MISet` containing the updated and sorted multi-indices.
"""
function add_mi(mi_set::MISet, mi_new::Vector)
    mi = copy(get_mi(mi_set))
    append!(mi, [mi_new])
    col_matrix = hcat(mi...)
    col_matrix_unique = unique(eachcol(col_matrix))
    col_matrix_sorted = sortmi(hcat(col_matrix_unique...))
    mi_set_new = MISet([Vector(v) for v in eachcol(col_matrix_sorted)])
    return mi_set_new
end

"""
    get_n_mi(mi_set)

Returns the number of multi-indices in a given `MISet`.

# Arguments
- `mi_set`: An instance of `MISet`.

# Returns
- The number of multi-indices in the set.
"""
function get_n_mi(mi_set::MISet)
    return length(mi_set.mi)
end

"""
    get_mi(mi_set)

Retrieves the list of multi-indices from a given `MISet`.

# Arguments
- `mi_set`: An instance of `MISet`.

# Returns
- A vector of multi-index vectors.
"""
function get_mi(mi_set::MISet)
    return mi_set.mi
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
    get_margin(mi_set)

Calculates the margin of a multi-index set (`mi_set`), which is the set of multi-indices that can extend the current set.

# Arguments
- `mi_set`: An instance of `MISet`.

# Returns
- An `MISet` containing the margin set.
"""
function get_margin(mi_set::MISet)
    mi = get_mi(mi_set)
    col_matrix = hcat(mi...)
    col_margin = margin(col_matrix)
    margin_set = MISet([Vector(v) for v in eachcol(col_margin)])
    return margin_set
end

"""
    get_reduced_margin(mi_set)

Calculates the reduced margin of a multi-index set (`mi_set`), which removes indices already in the set.

# Arguments
- `mi_set`: An instance of `MISet`.

# Returns
- An `MISet` containing the reduced margin set.
"""
function get_reduced_margin(mi_set::MISet)
    mi = copy(get_mi(mi_set))
    col_matrix = hcat(mi...)
    col_reduced_margin = reducedmargin(col_matrix)
    reduced_margin_set = MISet([Vector(v) for v in eachcol(col_reduced_margin)])
    return reduced_margin_set
end

"""
    create_smolyak_miset(n, k)

Creates a Smolyak multi-index set for a given dimension (`n`) and level (`k`).

# Arguments
- `n`: The dimensionality of the space.
- `k`: The level of the Smolyak grid.

# Returns
- An `MISet` representing the Smolyak multi-index set.
"""
function create_smolyak_miset(n, k)
    miset = createsmolyakmiset(n, k)
    miset_vec = [Vector(v) for v in eachcol(miset)]
    miset_smolyak = MISet(miset_vec)
    return miset_smolyak
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
    target_points_matrix = hcat(target_points...)'
    f_on_target_points = interpolateonsparsegrid(sg, f_on_grid, target_points_matrix)
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
    integral = integrateonsparsegrid(sg, f_on_grid, precompute; vecdims=nothing)
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
function integrate_L2_on_sparsegrid(sg, f_on_grid, precomputed_lagrange_integrals; product=dot, precomputed_pairwise_norms=nothing)
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
function precompute_pairwise_norms(f_on_grid, product)
    pairwisenorms = computepairwisenorms(f_on_grid, product)
    return pairwisenorms
end

end