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
    plot!(CCPoints(); npts=ii)
end
xlabel!("Clenshaw-Curtis Points")
```

Knot functions can be mapped to other domains.
Perhaps we wish to use Clenshaw-Curtis points on a parameter domain $\Gamma=[0,100]$.
```@example cc
a = 0
b = 100
p = plot(CCPoints(a,b))
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
    plot!(CCPoints(); npts=Doubling()(ii))
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
plot!(CCPoints(); npts=5)
plot!(UniformPoints(); npts=5)
plot!(gausschebyshevt(5)..., label="Gauss-Chebyshev t")
plot!(gausslegendre(5)..., label="Gauss-Legendre")
```
## Function Reference
```@autodocs
Modules = [SparseGridsKit]
Pages   = ["knots_structures.jl"]
```