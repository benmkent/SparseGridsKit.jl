```@example
using SparseGridsKit
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
                xlabel="y_1",
                ylabel="y_2",
                zlabel="y_3",
                camera = (20 * (1 + cos(i)),10 * (1 + cos(i)))
                )
end
```