# SparseGridsKit
A simple implementation of sparse grid polynomial interpolation and corresponding interpolation.

The construction is based upon the sparse operators introduced by Smolyak [Smolyak1963](@cite).
Key papers defining the model formulation include [Barthelmann2000](@cite), [Gerstner1998](@cite) and [Gerstner2003](@cite).
This construction of the approximation is inspired by the [Sparse Grids MATLAB Kit](https://sites.google.com/view/sparse-grids-kit) [Piazzola2024](@cite).

This documentation is example driven, with full documentation for each subset of functions given at the end of the relevant example.
This is not a thorough mathematical description of sparse grid approximation, but a practical guide to using the software.

The documentation contains:
```@contents
Depth = 1
```

The package is still under development and by no means complete.

An example sparse grid construction is illustrated below. 
```@example mixed
using SparseGridsKit, Plots, LaTeXStrings
# Test create_sparsegrid
n,k =3,3
knots = [linearpoints, ccpoints, uniformpoints]
rules = [linear, doubling, doubling]
mi_set = create_smolyak_miset(n,k)
sg = create_sparsegrid(mi_set, knots=knots, rule=rules)
nsteps = 100
@gif for i in range(0, stop = 2π, length = nsteps)
        plot_sparse_grid(sg)
        plot!(
                title="Sparse Grid n,k="*string(n)*","*string(k),
                xlabel=L"y_1",
                ylabel=L"y_2",
                zlabel=L"y_3",
                camera = (20 * (1 + cos(i)),10 * (1 + cos(i)))
                )
end
```