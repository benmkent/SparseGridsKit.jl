
# Adaptive Sparse Grids
Sparse grids can be adaptively constructed using a Gerstner--Griebel [Gerstner2003](@cite) type construction.
This is implemented in the [`adaptive_sparsegrid`](@ref) function.

## One dimensional example
To illustrate the adaptive algorithm consider functions that can be *exactly* approximated using a sparse grid approximation.
Consider the polynomial $f:[-1,1]\mapsto\R$
```math
f(x) = 3x^2 + 2x +1.
```
```@example 1d
    using SparseGridsKit
    nparams = 1
    f(x) = @. 3.0*x[1]^2 + 2.0*x[1] + 1.0
```
The [`adaptive_sparsegrid`](@ref) function can be called with a function $f$ and the dimension of the function domain `nparams`.
```@example 1d
(sg, f_on_Z) = adaptive_sparsegrid(f, nparams)
```
The function `f` can be exactly represented by a three point interpolant, which is identified in the adaptive algorithm.
```@example 1d
test_points = range(-1,stop=1,length=100)
all([f(x)] ≈ interpolate_on_sparsegrid(sg,f_on_Z,x) for x in test_points),  get_n_grid_points(sg) == 3
```
Taking powers $k$ of the polynomial $f$ gives a polynomial $f^k$ of polynomial degree $2^k$. Consequently, this is represented exactly using $2^k+1$ interpolation points.
The adaptive algorithm identifies this.
```@example 1d
        f2(x) = f(x).^2
        (sg, f2_on_Z) = adaptive_sparsegrid(f2, nparams)
        # Expect three point approx cubic (1 iteration to suffice)
        all([f2(x)] ≈ interpolate_on_sparsegrid(sg,f2_on_Z,x) for x in test_points), 
        get_n_grid_points(sg) == 5
```
```@example 1d
        f3(x) = f(x).^3
        (sg, f3_on_Z) = adaptive_sparsegrid(f3, nparams)
        # Expect three point approx cubic (1 iteration to suffice)
        all([f3(x)] ≈ interpolate_on_sparsegrid(sg,f3_on_Z,x) for x in test_points)
        get_n_grid_points(sg)
```
## Multi-dimensional example
The adaptive algorithm is also applicable for higher dimensional functions.
Consider the [`genz`](@ref) `gaussianpeak` example.
```@example genz
        using SparseGridsKit, LinearAlgebra
        n = 8
        C = 1.0
        W = 0.0
        T = "exponentialdecay"
        N = "gaussianpeak"
        f = genz(n::Int, C::Float64, W::Float64, T::String, N::String)
```
An adaptive approximation is formed.
This produces a (relatively) accurate approximation of the high-dimensional function $f$.
The iterations halt when the stopping criterion of max profit `p_max < proftol`.
```@example genz
    # Approximate
    (sg, f_on_Z) = adaptive_sparsegrid(f, n)

    nmc = 100
    test_points = [2 * rand(8) .- 1 for ii = 1:nmc]

    f_approx_on_test = interpolate_on_sparsegrid(sg,f_on_Z,test_points)

    norm([f(x) - f_approx_on_test[i] for (i,x) in enumerate(test_points)])
```
Reducing `proftol` allows further iterations and a more accurate approximation.
```@example genz
    (sg, f_on_Z) = adaptive_sparsegrid(f, n, proftol=1e-5)

    f_approx_on_test = interpolate_on_sparsegrid(sg,f_on_Z,test_points)

    norm([f(x) - f_approx_on_test[i] for (i,x) in enumerate(test_points)])
```

## Functions
```@autodocs
Modules = [SparseGridsKit]
Pages   = ["adaptivesparsegrids.jl"]
```