# Sparse Grids
```@example sg
using SparseGridsKit, Plots

n, k = 2, 1
mi_set = create_smolyak_miset(n, k)
sg = create_sparsegrid(mi_set)
propertynames(sg)
```

```@example sg
points = get_grid_points(sg)

x = [p[1] for p in points]
y = [p[2] for p in points]

scatter(x, y,
        title="Sparse Grid n,k="*string(n)*","*string(k),
        xlabel="y_1",
        ylabel="y_2")
```

The grid can be altered by adding additional multi-indices.
```@example sg
mi_set_new = MISet([[1,3]]) 
combined_mi_set = add_mi(mi_set, mi_set_new)
sg_combined = create_sparsegrid(combined_mi_set)

points = get_grid_points(sg_combined)

x = [p[1] for p in points]
y = [p[2] for p in points]

scatter(x, y,
        title="Sparse Grid Combined",
        xlabel="y_1",
        ylabel="y_2")
```
More complex grids can be constructed by increasing the dimension $n$ and adding multi-indices to the multi-index set. 
```@example sg
n, k = 3, 3
mi_set = create_smolyak_miset(n, k)
sg = create_sparsegrid(mi_set)
points = get_grid_points(sg)

x = [p[1] for p in points]
y = [p[2] for p in points]
z = [p[3] for p in points]

nsteps = 100
@gif for i in range(0, stop = 2π, length = nsteps)
        scatter(x, y, z,
                title="Sparse Grid n,k="*string(n)*","*string(k),
                xlabel="y_1",
                ylabel="y_2",
                zlabel="y_3",
                camera = (20 * (1 + cos(i)),10 * (1 + cos(i))))
end
```

## Function Reference
```@autodocs
Pages   = ["SparseGridsKit.jl"]
```
