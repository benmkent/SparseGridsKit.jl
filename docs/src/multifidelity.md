# Multi-fidelity Modelling
The construction offers support for multi-fidelity modelling.
A special type of knot function is included which simply maps a level to the same integer, to then be use to specify the model.
This has a weight $w=1.0$ and a corresponding polynomial equal to $1$ for all inputs.
In practice, for interpolation there is no parameter to be evaluated, but this set up provides the correct multi-fidelity structure using the existing single fidelity construction.
```@example mf
using SparseGridsKit
fidelitypoints(3), fidelity(3)
```

A simple multi-fidelity model considers a constant function $f(x)=y_1^2 + sin(y_2) + y_3$, subject to errors $10^-alpha$.
Functions are specified a $f(\vec{\alpha},\vec{y})$ where the vector $\vec{\alpha}$ controls the fidelity and $\vec{y}$ is simply the parameter as usual.
```@example mf
f(y) = y[1]^2 + sin(y[2]) + y[3]

nfid = 1
ndims = 3
f_fidelities(alpha,y) = f(y) + 10^(-alpha[1])
```

The domain is defined, where the fidelity domain should represent the maximum allowable fidelities.
The rule for fidelity is (`fidelity`)(@ref).
This returns $1$ for all input levels as there is only one model per level.
The knots are [`fidelitypoints`](@ref), which returns the input value plus a weight zero.
The knot function is treated as a special case in the sparse grid construction.
```@example mf
maxfidelity = 5
domain = [fill([1,maxfidelity],nfid)..., fill([-1,1],ndims)...]
rule = [fill(fidelity,nfid)..., fill(doubling,ndims)...]
knots = [fill(fidelitypoints,nfid)..., fill(ccpoints,ndims)...]
```

We wrap the function $f$ using [`multifidelityfunctionwrapper`](@ref).
This allows the user to use the single fidelity sparse grid construction, with function calls separated as $z=[\vec{\alpha},\vec{y}]$.
```@example mf
f_wrapped = multifidelityfunctionwrapper(f_fidelities,knots)

MI = create_smolyak_miset(nfid+ndims,2)

sg = create_sparsegrid(MI, domain; knots=knots, rule=rule)
@show get_grid_points(sg)
@show f_eval = f_wrapped.(get_grid_points(sg))
```
Finally, interpolation and integration are possible on the sparse grid.
Interpolation currently requires dummy values to give a complete parameter vector $[\vec{alpha} \vec{y}]$.
The values must be in the domain specified for the fidelity dimensions.
In the future this may be simplified.
```@example mf
interpolate_on_sparsegrid(sg, f_eval, [[fill(1,nfid)..., fill(0.0,ndims)...]])

pcl = precompute_lagrange_integrals(5, domain, knots, rule)

f_on_grid = f_wrapped.(get_grid_points(sg))

# Test integrate_on_sparsegrid
integral_result = integrate_on_sparsegrid(sg, f_on_grid, pcl)
```

It is also possible to use the adaptive sparse grid construction.
A `costfunction` can be supplied to give a representative cost for each fidelity.
This acts on a vector $\vec{α}$ of the fidelities.
```@example mf
(sg, f_on_Z) = adaptive_sparsegrid(f_wrapped, domain, ndims; maxpts = 100, proftol=1e-4, rule = rule, knots = knots, θ=1e-4, type=:deltaint, costfunction=α->10^prod(α))
```

## Function Reference
```@autodocs
Modules = [SparseGridsKit]
Pages   = ["misets.jl"]
```