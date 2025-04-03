# Plots

Plotting functionality is provided via [`Plots.jl](https://docs.juliaplots.org/latest/) and its `Recipies`.

## Knots
For the knot functions inheriting from the abstract type `Points` plotting is supported.
This defaults to 11 points.
```@example plots
using SparseGridsKit, Plots

p = plot(CCPoints(-3,3))
plot!(CCPoints(-3,3); npts=21)
```

## Multi-Index Set
For multi-index type `MISet` pairwise scatter plots for each dimension can be produced.
```@example plots
miset = create_smolyak_miset(3,3)
p = misetplot(miset)
```

## Sparse Grids
Sparse grid points in a `SparseGrid` can be plotted.
This has an optional arguement to select which indices are plotted which is useful in dimensions greater than 3.
```@example plots
sg = create_sparsegrid(miset,fill([-1,1],3))
p = plot(sg; targetdims=[3,2,1])
```

## Sparse Grid Approximations
Sparse grid approximations, consisting of function values and a polynomial approximation space can also be plotted.
```@example plots
f(x) = x[1]^5 + cos(x[2]) + abs(x[3])
f_on_sg = f.(get_grid_points(sg))

sga = SparseGridApproximation(sg,f_on_sg)
ssg = convert_to_spectral_approximation(sga)

p = plot(
    plot(sga; fill=true),
    plot(sga; targetdims=[2,3]),
    plot(ssg; color=:turbo, fill=true),
    plot(ssg; seriestype=:surface, targetdims=[2,3]),
    layout = (2,2)
)
```

## Function Reference
```@autodocs
Modules = [SparseGridsKit]
Pages   = ["sparsegridplots.jl"]
```