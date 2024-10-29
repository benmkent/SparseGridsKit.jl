# Sparse Grid Interpolation
## One-dimensional interpolation
Consider a simple example approximating a function $f:[-1,1]\to\R$.
```@example interp1
f(x) = [@. 3.0 * x[1]^2 + 2.0 * x[1] + 1.0]
```
We can construct a one point approximation (i.e. constant polynomial).
The sparse grid is constructed using the `MISet` `mi_set`.
The function `f` is evaluated on the grid points accessed via `get_grid_points`.
```@example interp1
n, k = 1, 0
using SparseGridsKit, LinearAlgebra
mi_set = create_smolyak_miset(n, k)
sg = create_sparsegrid(mi_set)
f_on_grid = [f(x) for x in get_grid_points(sg)]
```
The sparse grid `sg` is paired with the point evaluations to interpolate on to target points.
```@example interp1
target_points = [[x] for x in -1.0:0.5:1.0]
interpolation_result = interpolate_on_sparsegrid(sg, f_on_grid, target_points)
```
It is clear a better approximation can be formed.
```@example interp1
mi_set_new = create_smolyak_miset(n, 2)
sg_new = create_sparsegrid(mi_set_new)
f_on_grid_new = [f(x) for x in get_grid_points(sg_new)]

interpolation_result_new = interpolate_on_sparsegrid(sg_new, f_on_grid_new, target_points)

norm([interpolation_result_new[i] - f(target_points[i]) for i in eachindex(target_points)])
```
## Multi-dimensional interpolation
These ideas extend to multi-dimensional, vector valued functions.
Consider a function $f:[-1,1]^4 \to \mathbb{R}{400}$.
```@example interp4
using SparseGridsKit, LinearAlgebra
n, k = 4, 3
mi_set = create_smolyak_miset(n, k)
sg = create_sparsegrid(mi_set)
# Define a complicated function
f = genz(n, 1.0, 0.5, "quadraticdecay", "gaussianpeak")
f_on_grid = [[f(x)] for x in get_grid_points(sg)]
```
The sparse grid approximation defines a polynomial that can be evaluated.
```@example interp4
nmc = Integer(1e2);
V = [2 * rand(n) .- 1 for ii = 1:nmc]
f_on_v = [interpolate_on_sparsegrid(sg, f_on_grid, [v])[1] for v in V]
norm([f_on_v[i] - [f(v[i])] for i in eachindex(v)])
```
In this case this is poor.
We consider approximating the polynomial approximation.
If we grow the polynomial approximation space by adding multi-indices to the multi-index set we have a different polynomial approximation.
Consider the enriched approximation space formed by the union with the reduced margin.
```@example interp4
mi_enriched = add_mi(mi_set, get_reduced_margin(mi_set))
sg_enriched = create_sparsegrid(mi_enriched)
display(get_n_grid_points(sg_enriched))
```
Interpolating the polynomial via this new sparse grid should give the same function.
```@example interp4
f_on_grid_2 = interpolate_on_sparsegrid(sg, f_on_grid, get_grid_points(sg_enriched))
f2_on_v = interpolate_on_sparsegrid(sg_enriched, f_on_grid_2, v)
norm([f_on_v[i] - f2_on_v[i] for i in eachindex(v)])
```

The function $f$ is analytic so we expect the polynomial approximation to converge.
The approximation can be improved by increasing the size of the approximation space, at greater computational cost.
```@example interp4
n, k = 4, 6
mi_set = create_smolyak_miset(n, k)
sg = create_sparsegrid(mi_set)
f_on_grid = [f(x) for x in get_grid_points(sg)]
f_on_v = interpolate_on_sparsegrid(sg, f_on_grid, v)
norm([f_on_v[i] - [f(v[i])] for i in eachindex(v)])
```