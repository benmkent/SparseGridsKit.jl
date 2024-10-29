# Multi Index Sets

Multi-index sets the standard Smolyak construction are easily constructed.
Suppose we consider a function $u$ with domain $\Gamma\in\R^n$.
The multi-index set
$$ \{\alpha : \Vert \alpha \Vert_1 \leq n + d\} \subset \N^{2} $$
is constructed using the following syntax.
```@example misets
using SparseGridsKit

n, k = 2, 1
mi_set = create_smolyak_miset(n, k)
```
The multi-index set can be viewed a vector of `Vector{Integer}` representing each multi-index:
```@example misets
get_mi(mi_set)
```
The multi-index set can be altered by adding additional multi-indices.
```@example misets
mi_set_new = MISet([[1,3]]) 
combined_mi_set = add_mi(mi_set, mi_set_new)
```

In Gerstner--Griebel style algorithms, the margin and reduced margin are often required.
These are easily computed.
For example the margin,
```@example misets
margin_set = get_margin(combined_mi_set)
```
and the reduced margin,
```@example misets
reduced_margin_set = get_reduced_margin(combined_mi_set)
```
Notice how the multi-index $[2,3]$ is in the margin but not the reduced margin.
Its backwards neighbour $[2,2]$ is missing from the multi-index set.