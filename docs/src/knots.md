# Knots
## Knots functions
The knot functions are of the form
```
x,w = knots(n)
```
where `n` is the number of knots required.
The outputs are:
- `x` the vector of knots,
- and `w` the vector of quadrature weights.
The quadrature weights are generally normalised to sum to one (i.e. representing a probability measure).

Often, one chooses to use Clenshaw--Curtis points which are available using the function [`CCPoints()`](@ref).
```@example cc
using SparseGridsKit, Plots
p = plot()
for ii in 1:2:7
    x, w = CCPoints()(ii)
    scatter!(x, ii * ones(length(x)), label="n = $ii")
end
xlabel!("Clenshaw-Curtis Points")
```

Knot functions can be mapped to other domains.
Perhaps we wish to use Clenshaw-Curtis points on a parameter domain $\Gamma=[0,100]$.
```@example cc
a = 0
b = 100
x,w = CCPoints(a, b)(3)
p = plot()
scatter!(x, 0.0 * ones(length(x)))
xlabel!("Clenshaw-Curtis Points")
```

## Level functions
Additionally there are level functions mapping from a *level* $l$ to a number of knots $n$.
For example, we often use the approximate *doubling rule*
```math
n(l) = 2^{l-1} + 1 \quad \forall l > 1
```
with $n(1) = 1$.
```@example doubling
using SparseGridsKit
Doubling()(1), Doubling()(2), Doubling()(3)
```

Note that *sparsity* in the construction is only attained when using nested sets of points.
For Clenshaw-Curtis points this is achieved using the doubling rule.
```@example doubling
using SparseGridsKit, Plots
p = plot()
for ii in 1:5
    n = Doubling()(ii)
    x, w = CCPoints()(n)  # Use `ii` to vary the number of points
    scatter!(x, ii * ones(length(x)), label="n = $n")
end
xlabel!("Clenshaw-Curtis Points")
ylabel!("Level")
```

## More knot functions
Functionality is not restricted to knot functions included this package.
For example, one could use the quadrature points provided in the [`FastGaussQuadrature.jl`](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) package.
```@example
using SparseGridsKit, FastGaussQuadrature, Plots
p = plot()
x, w = CCPoints()(5)
scatter!(x, 1 * ones(length(x)), label="CC")
x, w = UniformPoints()(5)
scatter!(x, 2 * ones(length(x)), label="Uniform")
x, w = gausschebyshevt(5)
scatter!(x, 3 * ones(length(x)), label="Gauss-Chebyshev t")
x, w = gausslegendre(5)
scatter!(x, 4 * ones(length(x)), label="Gauss-Legendre")
```
## Function Reference
```@autodocs
Modules = [SparseGridsKit]
Pages   = ["knots_structures.jl"]
```