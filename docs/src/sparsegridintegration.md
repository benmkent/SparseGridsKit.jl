# Sparse Grid Integration
## Sparse Grid Quadrature
To integrate a sparse grid polynomial approximation we require the sparse grid, function evaluations and the qudarature weights.
These are precomputed when the sparse grid is constructed via the [`compute_quadrature_weights!`](@ref) function.

Consider the function $f:[-1,1]\to\mathbb{R}$
```@example int1
f(x) = @. 3.0*x^2 + 2.0*x + 1.0
```
and we construct the two-dimension function $[f(x_1), f(x_2)^2]^{\top}$ mapping $[-1,1]^2 \to \mathbb{R}^2$.
We seek the value of the integral
```math
\int_{\Gamma} [f(x_1), f(x_2)^2]^{\top} \rho(x) \textrm{d} x
```
for a weight function $\rho(x)=0.5^2$.

The constant approximation is $f(x)\approx[1,1]$.
The approximate integral using the constant approximation is equal to $[1,1]$.
```@example int1
using SparseGridsKit
n,k = 2,0
mi_set = create_smolyak_miset(n,k)

sg = create_sparsegrid(mi_set)
f_on_grid = [[f(x[1]), f(x[2])^2] for x in get_grid_points(sg)]
# Test integrate_on_sparsegrid
integral_result = integrate_on_sparsegrid(sg, f_on_grid)
result = [1, 1]
println("Quadrature: $integral_result, Integral: $result")
```
The estimated integral is improved by using a sparse grid with more grid points.
The integral can be evaluated explicitly, and it is seen that using an larger sparse grid gets the correct result $[2.0,92/15]$.
```@example int1
n,k = 2,4
mi_set = create_smolyak_miset(n,k)
sg = create_sparsegrid(mi_set)
f_on_grid = [[f(x[1]), f(x[2])^2] for x in get_grid_points(sg)]
# Test integrate_on_sparsegrid
integral_result = integrate_on_sparsegrid(sg, f_on_grid)
result = [2.0,92/15]
println("Quadrature: $integral_result, Integral: $result")
```
This is of course expected --- the polynomial approximation space represented by the larger sparse grid contains the function 
$[f(x_1), f(x_2)^2]$.

Underlying the integral is the computation of quadrature weights which are computed using [`compute_quadrature_weights!`](@ref).
This is done automatically when a grid is constructed.
```@example int1
compute_quadrature_weights!(sg)
sg.quadrature_weights
```

## Sparse Grid Weighted $L_{\rho}^2(\Gamma)$ Norms
Often the weighted $L_{\rho}^2(\Gamma)$ norm is required.
Squaring a polynomial approximation results in many higher degree polynomial terms --- it is not as simple as squaring the evaluations at each sparse grid point.

To approximate this integral this efficiently, integrals of the products of interpolation polynomials are precomputed.
A maximum level number is selected and for each domain dimension all pairwise $L_{\rho}^2(\Gamma)$ integrals inner products are computed using [`precompute_lagrange_integrals`](@ref).

For example, the pairwise integrals are computed, and we can extract the pairwise inner products for the level $2$ and level $3$ polynomials.
```@example int1
    using SparseGridsKit, LinearAlgebra
    pcl = precompute_lagrange_integrals(7)
    level1 = 2
    level2 = 3
    domaindim = 1
    pcl[domaindim][level1,level2,:,:]
```
This precomputation step is generally not too expensive but uses a moderate amount of memory.
```@example int1
@elapsed precompute_lagrange_integrals(7)
```
```@example int1
sizeof(pcl)
```
Additional arguments must be provided to  [`precompute_lagrange_integrals`](@ref) if specific knots and rules are use.

The precomputed $L_{\rho}^2(\Gamma)$ integrals from [`precompute_lagrange_integrals`](@ref) are used with the sparse grid approximation structure to *exactly* compute the $L_{\rho}^2(\Gamma)$ norm of a sparse grid polynomial function.
These are paired with appropriate pairwise products of the function evaluation which are computed with [`precompute_pairwise_norms`](@ref).
For example, for a scalar valued function we may require
```math
f(x_1) f(x_2) \forall x_1,x_2 \in Z
```
where $Z$ is the set of sparse grid points.

For example, the function $f(x)=x$ has the weighted $L_{\rho}^2([-1,1])$ norm equal to 
```math
\Vert x^2 \Vert_{L_{\rho}^2([-1,1])} = \sqrt{3}
```
for weight function $\rho=0.5$.
```@example int1
n,k = 1,1
mi_set = create_smolyak_miset(n,k)
pcl = precompute_lagrange_integrals(7)

sg = create_sparsegrid(mi_set)
f_on_grid = [x[1] for x in get_grid_points(sg)]
pairwise_norms = precompute_pairwise_norms(f_on_grid, product=(x,y)->dot(x,y))
l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
result = sqrt(1/3)
println("Quadrature: $l2_integral_result, Results: $result")
```
Similarly, for $x^2$ we get
```math
\Vert x^2 \Vert_{L_{\rho}^2([-1,1])}=\sqrt{5}.
``` 
```@example int1
n,k = 1,3
mi_set = create_smolyak_miset(n,k)
pcl = precompute_lagrange_integrals(7)

sg = create_sparsegrid(mi_set)
f_on_grid = [x[1]^2 for x in get_grid_points(sg)]
l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
result = sqrt(1/5)
println("Quadrature: $l2_integral_result, Results: $result")
```
Finally we consider the $L^2_{\rho}([-1,1]^2)$ norm of $f(x_1)$.
This can be explicitly computed to be $\sqrt{92/15}$.
```@example int1
n,k = 2,4
mi_set = create_smolyak_miset(n,k)
pcl = precompute_lagrange_integrals(7)

sg = create_sparsegrid(mi_set)
f_on_grid = [f(x[1]) for x in get_grid_points(sg)]
pairwise_norms = precompute_pairwise_norms(f_on_grid, product=(x,y)->x.*y)
l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
result = sqrt(92/15)
println("Quadrature: $l2_integral_result, Results: $result")
```

If the evaluation of the function $f:\Gamma\maps\to\R^n$ represents a discrete approximation of a function $F$, we may use a mass matrix $Q$ to compute
```math
(x_1^{\top} Q x_2)^{1/2} \quad \forall x_1,x_2 \in Z.
```
This is useful for integrals of the type $L^2_\rho(\Gamma; L^2(D))$ where $D$ represents the domain of the approximated function $F$ and is often seen in the approximation of parametric PDE solutions.
This is achieved by passing `product=(x,y)->dot(x,Q,y)` as an argument in [`precompute_pairwise_norms`](@ref). 

## Function Reference
```@docs
compute_quadrature_weights!
integrate_on_sparsegrid
integrate_L2_on_sparsegrid
```
