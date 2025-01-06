# Multi Index Sets
Sparse grids are constructed via *multi-index sets* specifying a cleverly chosen linear combination of tensor product grids.

Suppose we consider a function $u$ with domain $\Gamma\in\mathbb{R}^n$.

## Tensor Product Multi-Index Sets
The multi-index set
```math
\{\alpha : \Vert \alpha \Vert_\infty \leq n + d\} \subset \N^{d}
```
is constructed using the [`create_tensor_miset`](@ref) function.
```@example misets
using SparseGridsKit

n, k = 2, 1
mi_set = create_tensor_miset(n, k)
```
The multi-index set can be viewed a vector of `Vector{Integer}` representing each multi-index:
```@example misets
get_mi(mi_set)
```

## Smolyak Multi-Index Sets 
The standard Smolyak construction are also easily constructed.
The multi-index set
```math
\{\alpha : \Vert \alpha \Vert_1 \leq n + d\} \subset \N^{2}
```
is constructed using the following syntax.
```@example misets
using SparseGridsKit

n, k = 2, 1
mi_set = create_smolyak_miset(n, k)
```
This gives a different, smaller multi-index set than the tensor product multi-index set.
```@example misets
get_mi(mi_set)
```

## Manipulating Multi-Index Sets
The multi-index set can be altered by adding additional multi-indices.
Continuing using the above Smolyak example, we add the multi-index $[1, 3]$ with the [`add_mi`](@ref) function.
```@example misets
mi_set_new = MISet([[1,3]]) 
combined_mi_set = add_mi(mi_set, mi_set_new)
```

In Gerstner--Griebel style adaptive algorithms, the margin and reduced margin are generally required.
These are easily computed using the [`get_margin`](@ref) and [`get_reduced_margin`](@ref) functions.
```@example misets
margin_set = get_margin(combined_mi_set)
```
```@example misets
reduced_margin_set = get_reduced_margin(combined_mi_set)
```
Notice how the multi-index $[2,3]$ is in the margin but not the reduced margin.
Its backwards neighbour $[2,2]$ is missing from `combined_mi_set`.

## Function Reference
```@autodocs
Modules = [SparseGridsKit]
Pages   = ["misets.jl"]
```