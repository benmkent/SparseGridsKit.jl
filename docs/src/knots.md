# Knots
The knot functions take the form
```
x,w = knots(n)
```
where `n` is the number of knots required.
The outputs are:
- `x` the vector of knots,
- and `w` the vector of quadrature weights.
The quadrature weights are generally normalised to sum to one (i.e. representing a probability measure).

These underpin the sparse grid construction.

## Knot Functions Basics
Often, one chooses to use Clenshaw--Curtis points which are available using the function [`CCPoints()`](@ref).
```@example cc
using SparseGridsKit, Plots
p = plot()
for ii in 1:2:7
    plot!(CCPoints(); n=ii)
end
xlabel!("Clenshaw-Curtis Points")
```

Knot functions generally default to the unit interval $[-1,1]$.
They can also be mapped to other domains by passing an interval as an arguement.
For example, if we wish to use Clenshaw-Curtis points on a parameter domain $\Gamma=[0,100]$ we can do the following.
```@example cc
a = 0
b = 100
p = plot(CCPoints([a,b]))
xlabel!("Clenshaw-Curtis Points")
```

## Level Functions
Level functions allow us to define sequences of sets of knots with different growth rates.
The level functions are mappings from a *level* $l$ to a number of knots $n$.
For example, we often use the approximate *doubling rule*
```math
n(l) = 2^{l-1} + 1 \quad \forall l > 1, n(1) = 1.
```
```@example doubling
using SparseGridsKit
Doubling()(1), Doubling()(2), Doubling()(3)
```

Example level functions are plotted below.
```@example doubling
using Plots, LaTeXStrings
plot(1:5, Linear().(1:5),label="Linear")
plot(1:5, TwoStep().(1:5),label="TwoStep")
plot!(1:5, Doubling().(1:5),label="Doubling")
plot!(1:5, Tripling().(1:5),label="Tripling")
xlabel!(L"Level $l$")
ylabel!(L"$n(l)$")
```

Later, we will normally use nested sets of points to achieve *sparsity* in the sparse grid construction.
For Clenshaw-Curtis points, nestedness is achieved using the doubling rule.
```@example doubling
using SparseGridsKit, Plots
p = plot()
for ii in 1:5
    plot!(CCPoints(); n=Doubling()(ii))
end
xlabel!("Clenshaw-Curtis Points")
ylabel!("Level")
```

## More Knot Functions
Functionality is not restricted to knot functions included this package.
For example, one could use the quadrature points provided in the [`FastGaussQuadrature.jl`](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) package.
These can be wrapped as a [`CustomPoints`](@ref) to later be used to construct a sparse grid.
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

Leja points are also avaibile.
This currently either uses an discrete search based approach to iteratively construct points according to
```math
z^{k+1} = \textrm{argmax}_{z\in[-1,1]} v(z) \prod_{i=1}^{k} \textrm{abs}(z-z^i)
```
where $v(z)$ is the weight function.
The default weight is $v(z)= \sqrt(\rho(z))$ for $\rho(z)=0.5$ and nonsymmetric points.
```@example
using SparseGridsKit, Plots

p = plot(legend=:bottom)
symmetric = false
v(z) = sqrt(0.5)
plot!(LejaPoints(), label="default")
plot!(LejaPoints([-1,1],symmetric,:classic,v), label="nonsymmetric")
symmetric = true
plot!(LejaPoints([-1,1],symmetric,:classic), label="symmetric")
plot!(LejaPoints([-1,1],symmetric,:precomputed), label="precomputed,symmetric")
```

## Function Reference
```@autodocs
Modules = [SparseGridsKit]
Pages   = ["knots_structures.jl"]
```