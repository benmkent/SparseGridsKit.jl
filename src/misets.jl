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
    create_tensor_miset(n, k)

Creates a tensor product multi-index set for a given dimension (`n`) and level (`k`).

# Arguments
- `n`: The dimensionality of the space.
- `k`: The level of the tensor product grid.

# Returns
- An `MISet` representing the tensor product multi-index set.
"""
function create_tensor_miset(n, k)
    miset = createtensormiset(n, k)
    miset_vec = [Vector(v) for v in eachcol(miset)]
    miset_smolyak = MISet(miset_vec)
    return miset_smolyak
end