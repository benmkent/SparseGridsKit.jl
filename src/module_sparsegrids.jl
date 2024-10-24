module SparseGridsKit

export SparseGrid
export MISet

export create_sparsegrid, get_grid_points, get_n_grid_points, map_from_to, get_mi_set

export add_mi, get_mi, get_n_mi, get_margin, get_reduced_margin, create_smolyak_miset

export interpolate_on_sparsegrid, integrate_on_sparsegrid, integrate_L2_on_sparsegrid

export precompute_lagrange_integrals, precompute_pairwise_norms

## Use functionality defined in import file
include("sparsegrids.jl")

function precompute_lagrange_integrals(max_mi)
    precompute = sparsegridprecompute(max_mi)
    return precompute.productintegrals
end

function create_sparsegrid(mi_set)
    miset_matrix = hcat(mi_set.mi...)
    sg = createsparsegrid(miset_matrix; rule=doubling)
    return sg
end

function map_from_to(sg_from, sg_to)
    map_vector = mapfromto(sg_from, sg_to)
    return map_vector
end

function get_grid_points(sg)
    matrix_points = sg.sparsegridpts
    grid_point_vector = [Vector(v) for v in eachrow(matrix_points)]
    return grid_point_vector
end

function get_n_grid_points(sg)
    return size(sg.sparsegridpts, 1)
end

mutable struct MISet
    mi::Vector
end

function add_mi(mi_set::MISet, mi_set_new::MISet)
    mi = copy(get_mi(mi_set))
    mi_new = get_mi(mi_set_new)
    append!(mi, mi_new)
    col_matrix = hcat(mi...)
    col_matrix_sorted = sortmi(col_matrix)
    mi_set_new = MISet([Vector(v) for v in eachcol(col_matrix_sorted)])
    return mi_set_new
end

function add_mi(mi_set::MISet, mi_new::Vector)
    mi = copy(get_mi(mi_set))
    append!(mi, [mi_new])
    col_matrix = hcat(mi...)
    col_matrix_sorted = sortmi(col_matrix)
    mi_set_new = MISet([Vector(v) for v in eachcol(col_matrix_sorted)])
    return mi_set_new
end

function get_n_mi(mi_set::MISet)
    return length(mi_set.mi)
end

function get_mi(mi_set::MISet)
    return mi_set.mi
end

function get_mi_set(sg)
    mi = sg.MI
    mi_matrix = mi
    mi_matrix_dc = downwards_closed_set(mi_matrix)
    mi_set = MISet([Vector(v) for v in eachcol(mi_matrix_dc)])
    return mi_set
end

function get_margin(mi_set::MISet)
    mi = get_mi(mi_set)
    col_matrix = hcat(mi...)
    col_margin = margin(col_matrix)
    margin_set = MISet([Vector(v) for v in eachcol(col_margin)])
    return margin_set
end

function get_reduced_margin(mi_set::MISet)
    mi = copy(get_mi(mi_set))
    col_matrix = hcat(mi...)
    col_reduced_margin = reducedmargin(col_matrix)
    reduced_margin_set = MISet([Vector(v) for v in eachcol(col_reduced_margin)])
    return reduced_margin_set
end

function create_smolyak_miset(n, k)
    miset = createsmolyakmiset(n, k)
    miset_vec = [Vector(v) for v in eachcol(miset)]
    miset_smolyak = MISet(miset_vec)
    return miset_smolyak
end

function interpolate_on_sparsegrid(sg, f_on_grid, target_points)

    target_points_matrix = hcat(target_points...)'
    f_on_target_points = interpolateonsparsegrid(sg, f_on_grid, target_points_matrix)
    return f_on_target_points

end

function integrate_on_sparsegrid(sg, f_on_grid, precomputed_lagrange_integrals)
    precompute = (;productintegrals=precomputed_lagrange_integrals)

    integral = integrateonsparsegrid(sg, f_on_grid, precompute; vecdims=nothing)

    return integral
end

function integrate_L2_on_sparsegrid(sg, f_on_grid, precomputed_lagrange_integrals; product=dot, precomputed_pairwise_norms=nothing)
    precompute = (;productintegrals=precomputed_lagrange_integrals)
    integral = L2onsparsegrid(sg, f_on_grid, precompute; product=product, pairwisenorms=precomputed_pairwise_norms)
    return integral
end

function precompute_pairwise_norms(f_on_grid, product)
        pairwisenorms = computepairwisenorms(f_on_grid, product)
    return pairwisenorms
end
end