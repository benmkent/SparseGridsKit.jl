# Multi Index Sets
Sparse grids are constructed via *multi-index sets*.
These define a linear combination of tensor product grids, which are popular choices for approximation tasks in moderately high dimensional problems.

In the following, we consider a function $u$ with domain $\Gamma\in\mathbb{R}^n$.

## Tensor Product Multi-Index Sets
A tensor product multi-index set
```math
\{\alpha : \Vert \alpha \Vert_\infty \leq n + d\} \subset \N^{d}
```
is constructed using the [`create_tensor_miset`](@ref) function.
```@example misets
using SparseGridsKit, Plots, LaTeXStrings

n, k = 2, 1
mi_set = create_tensor_miset(n, k)
```
The multi-index set can be viewed a vector of `Vector{Integer}` representing each multi-index:
```@example misets
get_mi(mi_set)
```
```@example misets
x = getindex.(get_mi(mi_set), 1)
y = getindex.(get_mi(mi_set), 2)
scatter(x,y, xlabel=L"\beta_1", ylabel=L"\beta_2", 
            legend=:none,
            title="Multi-Index Set", aspect_ratio=:equal)
```
The exponential growth in complexity with domain dimension $n$ is problematic as soon as $n$ becomes moderately large.

## Smolyak Multi-Index Sets 
To alleviate some of the problems of high-dimensional approximation, Smolyak type constructions are often used.
The standard Smolyak multi-index sets are can be easily constructed.
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
```@example misets
x = getindex.(get_mi(mi_set), 1)
y = getindex.(get_mi(mi_set), 2)
scatter(x,y, xlabel=L"\beta_1", ylabel=L"\beta_2", 
            legend=:none,
            title="Smolyak Multi-Index Set", aspect_ratio=:equal)
```

## Manipulating Multi-Index Sets
The multi-index set can be altered by adding additional multi-indices.
Continuing using the above Smolyak example, we add the multi-index $[1, 3]$ with the [`add_mi`](@ref) function.
```@example misets
mi_set_new = MISet([[1,3]]) 
combined_mi_set = add_mi(mi_set, mi_set_new)

x = getindex.(get_mi(combined_mi_set), 1)
y = getindex.(get_mi(combined_mi_set), 2)
scatter(x,y, xlabel=L"\beta_1", ylabel=L"\beta_2", 
            legend=:none,
            title="Multi-Index Addition", aspect_ratio=:equal)
```

In Gerstner--Griebel style adaptive algorithms, the margin and reduced margin are generally required.
These are easily computed using the [`get_margin`](@ref) and [`get_reduced_margin`](@ref) functions.
```@example misets
margin_set = get_margin(combined_mi_set)

x = getindex.(get_mi(margin_set), 1)
y = getindex.(get_mi(margin_set), 2)
scatter(x,y, xlabel=L"\beta_1", ylabel=L"\beta_2", 
            legend=:none,
            title="Margin Multi-Index", aspect_ratio=:equal)
```
```@example misets
reduced_margin_set = get_reduced_margin(combined_mi_set)

x = getindex.(get_mi(reduced_margin_set), 1)
y = getindex.(get_mi(reduced_margin_set), 2)
scatter(x,y, xlabel=L"\beta_1", ylabel=L"\beta_2", 
            legend=:none,
            title="Reduced Margin Multi-Index", aspect_ratio=:equal)
```
Notice how the multi-index $[2,3]$ is in the margin but not the reduced margin.
Its backwards neighbour $[2,2]$ is missing from `combined_mi_set`.

## Function Reference
```@autodocs
Modules = [SparseGridsKit]
Pages   = ["misets.jl"]
```