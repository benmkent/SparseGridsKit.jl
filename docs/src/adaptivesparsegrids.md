
# Adaptive Sparse Grids
Sparse grids can be adaptively constructed using a Gerstner--Griebel [Gerstner2003](@cite) type construction.

## One dimensional example
To illustrate the adaptive algorithm consider functions that can be *exactly* approximated using the sparse grid approximation.

```@example 1d
    using SparseGridsKit
    f(x) = @. 3.0*x[1]^2 + 2.0*x[1] + 1.0
    ndims = 1
```
The `adaptive_sparsegrid` function can be called with a function $f$ and the dimension of the function domain `ndims`.
```@example 1d
(sg, f_on_Z) = adaptive_sparsegrid(f, ndims)
```
This can be exactly represented by a three point interpolant, which is identified in the adaptive algorithm.
```@example 1d
test_points = range(-1,stop=1,length=100)
display(all([f(x)] ≈ interpolate_on_sparsegrid(sg,f_on_Z,x) for x in test_points))
display(get_n_grid_points(sg) == 3)
```
Taking powers $k$ of the polynomial $f$ gives a polynomial $f^k$ of polynomial degree $2^k$. Consequently, this is approximated exactly using $2^k+1$ interpolation points.
The adaptive algorithm identifies this.
```@example 1d
        f2(x) = f(x).^2
        (sg, f2_on_Z) = adaptive_sparsegrid(f2, ndims)
        # Expect three point approx cubic (1 iteration to suffice)
        display(all([f2(x)] ≈ interpolate_on_sparsegrid(sg,f2_on_Z,x) for x in test_points))
        display(get_n_grid_points(sg) == 5)
```
```@example 1d
        f3(x) = f(x).^3
        (sg, f3_on_Z) = adaptive_sparsegrid(f3, ndims)
        # Expect three point approx cubic (1 iteration to suffice)
        display(all([f3(x)] ≈ interpolate_on_sparsegrid(sg,f3_on_Z,x) for x in test_points))
        display(get_n_grid_points(sg) == 9)
```
## Multi-dimensional example
The adaptive algorithm is also applicable for higher dimensional functions.
Consider the Genz `gaussianpeak` example.
```@example genz
        n = 8
        C = 1.0
        W = 0.0
        T = "exponentialdecay"
        N = "gaussianpeak"
        f = genz(n::Int, C::Float64, W::Float64, T::String, N::String)
```
An adaptive approximation is created.
This produces a (relatively) accurate approximation of the high-dimensional function $f$.
```@example
        # Approximate
    (sg, f_on_Z) = adaptive_sparsegrid(f, n)

    sg_test = create_sparsegrid(create_smolyak_miset(n,3))
    test_points = get_grid_points(sg_test)

    f_approx_on_test = interpolate_on_sparsegrid(sg,f_on_Z,test_points)

    display(all(isapprox(f(x), f_approx_on_test[i]; atol=1e-3) for (i,x) in enumerate(test_points)))
```

## Functions
```@autodocs
Pages   = ["adaptivesparsegrids.jl"]
```