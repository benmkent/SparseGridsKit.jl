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
    plot!(CCPoints(); n=ii)
end
xlabel!("Clenshaw-Curtis Points")
```

Knot functions can be mapped to other domains.
Perhaps we wish to use Clenshaw-Curtis points on a parameter domain $\Gamma=[0,100]$.
```@example cc
a = 0
b = 100
p = plot(CCPoints([a,b]))
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
    plot!(CCPoints(); n=Doubling()(ii))
end
xlabel!("Clenshaw-Curtis Points")
ylabel!("Level")
```

## More knot functions
Functionality is not restricted to knot functions included this package.
For example, one could use the quadrature points provided in the [`FastGaussQuadrature.jl`](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) package.
These can be wrapped in [`CustomPoints`](@ref).
```@example
using SparseGridsKit, FastGaussQuadrature, Plots

CustomGaussChebyshevPoints = CustomPoints([-1.0,1.0], n->gausschebyshevt(n))
CustomGaussLegendrePoints = CustomPoints([-1.0,1.0], n->gausslegendre(n))

p = plot()
plot!(CCPoints(); n=5)
plot!(UniformPoints(); n=5)
plot!(CustomGaussChebyshevPoints; n=5)
plot!(CustomGaussLegendrePoints; n=5)
```
Leja points are also supplied.
This currently either uses an discrete search based approach to iteratively construct points according to
```math
z^{k+1} = \argmax_{z\in[-1,1]} v(z) \prod_{i=1}^{k} \abs(z-z^i)
```
where $v(z)$ is the weight function, or a precomputed set for a uniform weight function.
<!-- Optimisation based points differ slightly to points generated using a discrete search. -->
The default weight is $v(z)= \sqrt(\rho(z))$ for $\rho(z)=0.5$ and unsymmetrical points.
```@example
using SparseGridsKit, Plots

p = plot()
symmetric = false
v(z) = sqrt(0.5)
plot!(LejaPoints(), label="default")
plot!(LejaPoints([-1,1],symmetric,:classic,v), label="unsymmetrical")
symmetric = true
plot!(LejaPoints([-1,1],symmetric,:classic), label="symmetrical")
plot!(LejaPoints([-1,1],symmetric,:precomputed), label="precomputed")
```

## Function Reference
```@autodocs
Modules = [SparseGridsKit]
Pages   = ["knots_structures.jl"]
```