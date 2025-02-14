# Sparse Grid Integration
Sparse grid interpolation in `SparseGridsKit.jl` relies on precomputed integrals of the underlying approximation functions.

!!! note
    This may change in the future - at the moment precomputing the pairwise integrals is used to speed up computation of $L^2(\Gamma)$ type norms.

To do this a maximum level number is selected and for each domain dimension, all pairwise $L_{\rho}^2(\Gamma)$ integrals inner products are computed using [`precompute_lagrange_integrals`](@ref).
For example, we can then extract the pairwise inner products for the level $2$ and level $3$ polynomials.
```@example int1
    using SparseGridsKit, LinearAlgebra
    domain = [[-1,1]]
    pcl = precompute_lagrange_integrals(7,domain)
    level1 = 2
    level2 = 3
    domaindim = 1
    pcl[domaindim][level1,level2,:,:]
```
This precomputation step is generally not too expensive but uses a moderate amount of memory.
```@example int1
@elapsed precompute_lagrange_integrals(7,domain)
```
```@example int1
sizeof(pcl)
```

Additional arguments must be provided to  [`precompute_lagrange_integrals`](@ref) if specific knot types are use.

## Sparse Grid Integration
To integrate a sparse grid polynomial approximation we require the sparse grid, function evaluations and the precomputed integrals.

!!! note
    The use of the full precomputed pairwise integrals is overkill here. We only require the weights output by each knot function. This should be optimised in the future.

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
    n,k = 2,0
    mi_set = create_smolyak_miset(n,k)
    domain = fill([-1,1],2)
    pcl = precompute_lagrange_integrals(7,domain)

    sg = create_sparsegrid(mi_set, domain)
    f_on_grid = [[f(x[1]), f(x[2])^2] for x in get_grid_points(sg)]
    # Test integrate_on_sparsegrid
    integral_result = integrate_on_sparsegrid(sg, f_on_grid, pcl)
```
The estimated integral is improved by using a sparse grid with more grid points.
The integral can be evaluated explicitly, and it is seen that using an larger sparse grid gets the correct result $[2.0,92/15]$.
```@example int1
    n,k = 2,4
    mi_set = create_smolyak_miset(n,k)
    domain = [[-1,1],[-1,1]]
    sg = create_sparsegrid(mi_set,domain)
    f_on_grid = [[f(x[1]), f(x[2])^2] for x in get_grid_points(sg)]
    # Test integrate_on_sparsegrid
    integral_result = integrate_on_sparsegrid(sg, f_on_grid, pcl)
```
This is of course expected --- the polynomial approximation space represented by the larger sparse grid contains the function 
$[f(x_1), f(x_2)^2]$.

## Sparse Grid Weighted L_{\rho}^2(\Gamma) Norms
Often the weighted $L_{\rho}^2(\Gamma)$ norm is required.
Squaring a polynomial approximation results in many higher degree polynomial terms --- it is not as simple as squaring the evaluations at each sparse grid point.

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
    domain = [[-1,1]]
    pcl = precompute_lagrange_integrals(7,domain)

    sg = create_sparsegrid(mi_set,domain)
    f_on_grid = [[x[1]] for x in get_grid_points(sg)]
    pairwise_norms = precompute_pairwise_norms(f_on_grid, product=(x,y)->dot(x,y))
    l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
    l2_integral_result[1] ≈ sqrt(1/3)
```
Similarly, for $x^2$ we get
```math
\Vert x^2 \Vert_{L_{\rho}^2([-1,1])}=\sqrt{5}.
``` 
```@example int1
    n,k = 1,3
    mi_set = create_smolyak_miset(n,k)
    domain = [[-1,1]]
    pcl = precompute_lagrange_integrals(7,domain)

    sg = create_sparsegrid(mi_set,domain)
    f_on_grid = [[x[1]^2] for x in get_grid_points(sg)]
    l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
    l2_integral_result[1] ≈ sqrt(1/5)
```
Finally we consider the $L^2_{\rho}([-1,1]^2)$ norm of $f(x_1)$.
This can be explicitly computed to be $\sqrt{92/15}$.
```@example int1
    n,k = 2,4
    mi_set = create_smolyak_miset(n,k)
    domain = [[-1,1],[-1,1]]
    sg = create_sparsegrid(mi_set,domain)
    f_on_grid = [f(x[1]) for x in get_grid_points(sg)]
    pairwise_norms = precompute_pairwise_norms(f_on_grid, product=(x,y)->x.*y)
    l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
    l2_integral_result ≈ sqrt(92/15)
```

If the evaluation of the function $f:\Gamma\maps\to\R^n$ represents a discrete approximation of a function $F$, we may use a mass matrix $Q$ to compute
```math
(x_1^{\top} Q x_2)^{1/2} \quad \forall x_1,x_2 \in Z.
```
This is useful for integrals of the type $L^2_\rho(\Gamma; L^2(D))$ where $D$ represents the domain of the approximated function $F$ and is often seen in the approximation of parametric PDE solutions.
This is achieved by passing `product=(x,y)->dot(x,Q,y)` as an argument in [`precompute_pairwise_norms`](@ref). 

## Function Reference
```@docs
integrate_on_sparsegrid
integrate_L2_on_sparsegrid
```
