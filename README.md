# SparseGridsKit.jl
[![Documentation](https://github.com/benmkent/SparseGridsKit.jl/actions/workflows/documentation.yaml/badge.svg)](http://benmkent.github.io/SparseGridsKit.jl/)
[![Test](https://github.com/benmkent/SparseGridsKit.jl/actions/workflows/test.yaml/badge.svg)](https://github.com/benmkent/SparseGridsKit.jl/actions/workflows/test.yaml)
[![codecov](https://codecov.io/github/benmkent/SparseGridsKit.jl/graph/badge.svg?token=URGWM64U21)](https://codecov.io/github/benmkent/SparseGridsKit.jl)

A simple implementation of sparse grid polynomial interpolation and corresponding interpolation.

The construction is based upon the sparse operators introduced by Smolyak [Smolyak1963](@cite).
Key papers defining the model formulation include [Barthelmann2000](@cite), [Gerstner1998](@cite) and [Gerstner2003](@cite).
This construction of the approximation is inspired by the [Sparse Grids MATLAB Kit](https://sites.google.com/view/sparse-grids-kit) [Piazzola2024](@cite).

The package is still under development and by no means complete.

# Example
In the simplest case one often wants to create a sparse grid using a sequence of nested sets Clenshaw--Curtis points and a Smolyak type multi-index set.
To this end, one uses `create_sparsegrid` with the desired Smolyak multi-index set.
In this case we choose a grid with dimension `n=3` and Smolyak level `w=3`.
```
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
![Example Sparse Grid](ExampleGrid.gif)
