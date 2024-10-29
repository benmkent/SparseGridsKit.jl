# Sparse Grid Integration
Sparse grid interpolation in `SparseGridsKit.jl` relies on precomputed integrals of the Lagrange interpolation polynomials.

Integrals are currently computed solely assuming Clenshaw--Curtis points and a density function $\rho=\frac{1}{2}^{n}$ where $n$ is the dimension of the parameter space $\Gamma=[-1,1]^{n}$

To do this a maximum level number is selected.
```@example int1
    using SparseGridsKit, LinearAlgebra
    pcl = precompute_lagrange_integrals(7)
    pcl[2,1,:,:]
```
## Sparse Grid Integration
This is then used to integrate a sparse grid polynomial approximation.
Consider the function $f:[-1,1]\to\R$
```@example int1
    f(x) = @. 3.0*x^2 + 2.0*x + 1.0
```
The two-dimension function $[f(x_1) f(x_2)^2]$ mapping $[-1,1]^2 \to \R^2$ can be approximated using a constant function.
The integral is expected (albeit a poor approximation to the true value).
```@example int1
    n,k = 2,0
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    f_on_grid = [[f(x[1]), f(x[2])^2] for x in get_grid_points(sg)]
    # Test integrate_on_sparsegrid
    integral_result = integrate_on_sparsegrid(sg, f_on_grid, pcl)
    integral_result ≈ [1.0,1.0]
```
The integral estimate is improved by using a sparse grid with more grid points.
The integral can be evaluated explicitly, and it is seen that the quadrature gets the correct result.
```@example int1
    n,k = 2,4
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    f_on_grid = [[f(x[1]), f(x[2])^2] for x in get_grid_points(sg)]
    # Test integrate_on_sparsegrid
    integral_result = integrate_on_sparsegrid(sg, f_on_grid, pcl)

    integral_result ≈ [2.0,92/15]
```

## Sparse Grid L^2(Gamma) Norms
Often the weighted $L^2(Gamma)$ norm is required.
Squaring a polynomial approximation results in many higher degree polynomial terms --- it is not as simple as squaring the evaluations at each sparse grid point.

The precomputed $L_{\rho}^2(Gamma)$ integrals in `pcl` are used with the sparse grid approximation structure to *exactly* compute the $L_{\rho}^2(\Gamma)$ norm of a sparse grid polynomial function.
This requires the pairwise norms of the point evaluations to be computed using `precompute_pairwise_norms`.
If the evaluations are vectors, a mass matrix can be supplied in `(x,y)->dot(x,Q,y)`. 

For the function `x` the weighted $L_{\rho}^2([-1,1])$ norm is $ $\Vert x^2 \Vert_{L_{\rho}^2([-1,1])} = \sqrt{3}'.
```@example int1
    n,k = 1,1
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    f_on_grid = [[x[1]] for x in get_grid_points(sg)]
    pairwise_norms = precompute_pairwise_norms(f_on_grid, (x,y)->dot(x,y))
    l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
    l2_integral_result[1] ≈ sqrt(1/3)
```
Similarly, for $x^2$ we get $\Vert x^2 \Vert_{L_{\rho}^2([-1,1])}=\sqrt{5}$. 
```@example int1
    n,k = 1,3
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    f_on_grid = [[x[1]^2] for x in get_grid_points(sg)]
    l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
    l2_integral_result[1] ≈ sqrt(1/5)
```
Finally we consider the $L^2_{\rho}([-1,1]^2)$ norm of $f(x_1)$.
This can be explicitly computed to be $\sqrt(92/15)$.
```@example int1
    # Test precompute_pairwise_norms
    n,k = 2,4
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set)
    f_on_grid = [[f(x[1])] for x in get_grid_points(sg)]
    pairwise_norms = precompute_pairwise_norms(f_on_grid, (x,y)->x.*y)
    l2_integral_result = integrate_L2_on_sparsegrid(sg, f_on_grid, pcl)
    l2_integral_result ≈ sqrt(92/15)
```