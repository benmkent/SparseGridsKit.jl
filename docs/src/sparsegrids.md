# Sparse Grids
Knot functions, level functions and multi-index sets are combined to define a sparse grid structure.

## Simple Sparse Grids
In the simplest case one often wants to create a sparse grid using a sequence of nested sets Clenshaw--Curtis points and a Smolyak type multi-index set.
To this end, one uses [`create_sparsegrid`](@ref) with the desired Smolyak multi-index set.
```@example sg
using SparseGridsKit

n, k = 2, 1
mi_set = create_smolyak_miset(n, k)
domain = [[-1,1],[-1,1]]
sg = create_sparsegrid(mi_set,domain)
propertynames(sg)
```
The grid points can be extracted using the [`get_grid_points`](@ref) function.
This returns a vector of vector points.
```@example sg
points = get_grid_points(sg)
```

```@example sg
using Plots, LaTeXStrings
plot_sparse_grid(sg)
plot!(title="Sparse Grid n,k="*string(n)*","*string(k),
        xlabel=L"y_1",
        ylabel=L"y_2")
```

The grid can be altered by adding additional multi-indices.
```@example sg
mi_set_new = MISet([[1,3]]) 
combined_mi_set = add_mi(mi_set, mi_set_new)
sg_combined = create_sparsegrid(combined_mi_set,domain)

points = get_grid_points(sg_combined)

x = [p[1] for p in points]
y = [p[2] for p in points]

plot_sparse_grid(sg_combined)
plot!(  title="Sparse Grid Combined",
        xlabel="y_1",
        ylabel="y_2")
```
More complex grids can be constructed.
For example consider increasing the dimension $n$ and using more multi-indices in the multi-index set. 
```@example sg
n, k = 3, 3
mi_set = create_smolyak_miset(n, k)
domain = [[-1,1],[-1,1],[-1,1]]
sg = create_sparsegrid(mi_set,domain)
points = get_grid_points(sg)
```
This can still be easily visualised.
```@example sg
nsteps = 100
@gif for i in range(0, stop = 2π, length = nsteps)
        plot_sparse_grid(sg)
        plot!(title="Sparse Grid n,k="*string(n)*","*string(k),
                xlabel=L"y_1",
                ylabel=L"y_2",
                zlabel=L"y_3",
                camera = (20 * (1 + cos(i)),10 * (1 + cos(i))))
end
```
## Mixed Knots Sparse Grids
It is possible to mix the type of points in a grid.
To demonstrate, consider linear points $$n\mapsto\{1,2,...,n\}$$ in the first dimension, Clenshaw--Curtis points on a domain $\Gamma_2=[0,100]$ in the second dimension and uniformly spaced points on the domain $[-1,1]$  in the third dimension.
```@example mixed
using SparseGridsKit, Plots, LaTeXStrings
# Test create_sparsegrid
n,k =3,3
knots = [linearpoints, n->ccpoints(n,0,100), uniformpoints]
rules = [linear, doubling, doubling]
mi_set = create_smolyak_miset(n,k)
domain = [[-1,1],[0,100],[0,1]]
sg = create_sparsegrid(mi_set, domain, knots=knots, rule=rules)
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

## Function Reference
```@autodocs
Modules = [SparseGridsKit]
Pages   = ["SparseGridsKit.jl"]
```
