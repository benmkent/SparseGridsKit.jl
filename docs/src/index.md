# SparseGridsKit.jl
A simple implementation of sparse grid polynomial interpolation and corresponding interpolation.
Its aim is to provide adaptive approximation code closely resembling the mathematical literature to aid algorithm development and analysis.

To add the package use the following command:
```
] add https://github.com/benmkent/SparseGridsKit.jl
```
The construction is based upon the sparse operators introduced by Smolyak [Smolyak1963](@cite).
Key papers defining the model formulation include [Barthelmann2000](@cite), [Gerstner1998](@cite) and [Gerstner2003](@cite).
This construction of the approximation is inspired by the [Sparse Grids MATLAB Kit](https://sites.google.com/view/sparse-grids-kit) [Piazzola2024](@cite).

This documentation is example driven, with full documentation for each subset of functions given at the end of the relevant example.
This is not a thorough mathematical description of sparse grid approximation, but a practical guide to using the software.

The package is still under development and by no means complete.

## Contributing
To report bugs, seek support or submit request features, please use the [GitHub issue tracker](https://github.com/benmkent/SparseGridsKit.jl/issues).
Please feel free to contribute and submit pull requests to the main branch.

## Example
Please see the [documentation](https://benmkent.github.io/SparseGridsKit.jl/stable/) for extensive examples.

An example sparse grid construction is illustrated below.
```@example mixed
using SparseGridsKit, Plots, LaTeXStrings
# Test create_sparsegrid
n,k =3,3
knots = [CCPoints(), CCPoints(), UniformPoints()]
rules = [Doubling(), Doubling(), Linear()]
mi_set = create_smolyak_miset(n,k)
sg = create_sparsegrid(mi_set, knots=knots, rule=rules)
nsteps = 100
@gif for i in range(0, stop = 2π, length = nsteps)
        plot(sg)
        plot!(
                title="Sparse Grid n,k="*string(n)*","*string(k),
                xlabel=L"y_1",
                ylabel=L"y_2",
                zlabel=L"y_3",
                camera = (20 * (1 + cos(i)),10 * (1 + cos(i)))
                )
end
```
Paired with function evaluations `f_on_Z = [f(z) for z in get_grid_points(sg)]` at the sparse grid points `Z = get_grid_points(sg)`, this defines a sparse grid approximation.
