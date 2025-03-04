# Sparse Grid Interpolation
Once a sparse grid is constructed, a polynomial approximation is formed by pairing the grid points with function evaluations.

If the multi-index set is admissible, and the sequences of sets of points are nested this approximation is an interpolating polynomial.

## One-dimensional interpolation
Consider a simple example approximating the function $f:[-1,1]\to\mathbb{R}$
```math
f(x) = 3x^2 + 2x + 1
```
Ths is expressed as the function `f`.
```@example interp1
f(x) = [@. 3.0 * x[1]^2 + 2.0 * x[1] + 1.0]
```
Notice that the function `f` returns a one element `Vector`.
This is not essential, but is advisable as vectors and scalars are treated differently in `julia`.

We can construct a one point approximation (i.e. constant polynomial).
The function `f` is evaluated on the grid points accessed via [`get_grid_points`](@ref).
```@example interp1
n, k = 1, 0
using SparseGridsKit, LinearAlgebra
mi_set = create_smolyak_miset(n, k)
domain = [[-1,1]]
sg = create_sparsegrid(mi_set,domain)
f_on_grid = [f(x) for x in get_grid_points(sg)]
```
The sparse grid `sg` is paired with the point evaluations to interpolate on to target points.
Notice that whilst we are in one-dimension, `target_points` is still a `Vector` of `Vector`s
```@example interp1
target_points = [[x] for x in -1.0:0.5:1.0]
interpolation_result = interpolate_on_sparsegrid(sg, f_on_grid, target_points)
```
The constant approximation is $f\approx 1$.
It is clear a better approximation can be formed.
This is attained with a higher degree polynomial approximation -- derived from a larger multi-index set.
```@example interp1
mi_set_new = create_smolyak_miset(n, 1)
sg_new = create_sparsegrid(mi_set_new,domain)
f_on_grid_new = [f(x) for x in get_grid_points(sg_new)]

interpolation_result_new = interpolate_on_sparsegrid(sg_new, f_on_grid_new, target_points)

norm([interpolation_result_new[i] - f(target_points[i]) for i in eachindex(target_points)]) < 1e-5
```

The number of grid points can be checked using [`get_n_grid_points`](@ref) and is seen to equal three.
It is therefore no surprise that the approximation error is at machine precision - we can exactly represent the function `f` as a three point interpolant.
```@example interp1
get_n_grid_points(sg_new) == 3
```
## Multi-dimensional interpolation
These ideas extend to vector valued functions with multi-dimensional domains.
Consider a function
```math
f:[-1,1]^4 \to \mathbb{R}^{400}
```
defined via the `gaussianpeak` Genz test function (see e.g. [Barthelmann2000](@cite)).
A polynomial approximation can again be formed using nested Clenshaw-Curtis points and Smolyak multi-index sets.

```@example interp4
using SparseGridsKit, LinearAlgebra, Plots, LaTeXStrings
n, k = 4, 3
mi_set = create_smolyak_miset(n, k)
domain = [[-1,1],[-1,1],[-1,1],[-1,1]]
sg = create_sparsegrid(mi_set,domain)

# Define a complicated function
ndims=100
f(x) = [genz(n, 1.0, w, "quadraticdecay", "gaussianpeak")(x) for (i,w) = enumerate(range(-1,stop=1,length=ndims))]
f_on_grid = [f(x) for x in get_grid_points(sg)]

xplot = range(-1,stop=1,length=100)
plot(xplot, [f([xi, zeros(3)...])[1] for xi in xplot],label=L"f_{1}(x)")
plot!(xplot, [f([xi, zeros(3)...])[25] for xi in xplot],label=L"f_{25}(x)")
plot!(xplot, [f([xi, zeros(3)...])[50] for xi in xplot],label=L"f_{50}(x)")
plot!(xplot, [f([xi, zeros(3)...])[100] for xi in xplot],label=L"f_{100}(x)")
plot!(
    title=L"f(x_1,0,0,...)",
    xlabel=L"x_1"
)
```
This sparse grid approximation is constructed and evaluated at $100$ test points.
The norm of the error vector is computed and is small, but certainly not at machine precision.
```@example interp4
nmc = Integer(1e2);
V = [2 * rand(n) .- 1 for ii = 1:nmc]
f_on_v = [interpolate_on_sparsegrid(sg, f_on_grid, [v])[1] for v in V]
norm([f_on_v[i] - f(V[i]) for i in eachindex(V)])
```
If we add multi-indices to the multi-index set we should have a better, higher degree polynomial approximation.
The function $f$ is analytic so we expect the polynomial approximation to converge quickly.
```@example interp4
n, k = 4, 6
mi_set = create_smolyak_miset(n, k)
sg = create_sparsegrid(mi_set,domain)
f_on_grid = [f(x) for x in get_grid_points(sg)]
f_on_v = interpolate_on_sparsegrid(sg, f_on_grid, V)
norm([f_on_v[i] - f(V[i]) for i in eachindex(V)])
```
```@example interp4
n, k = 4, 7
mi_set = create_smolyak_miset(n, k)
sg = create_sparsegrid(mi_set,domain)
f_on_grid = [f(x) for x in get_grid_points(sg)]
f_on_v = interpolate_on_sparsegrid(sg, f_on_grid, V)
norm([f_on_v[i] - f(V[i]) for i in eachindex(V)]) < 1e-5
```
Later, the function [`adaptive_sparsegrid`](@ref) will consider a greedy construction of the multi-index set to reduce the approximation error.

## Function Reference
```@docs
interpolate_on_sparsegrid
```