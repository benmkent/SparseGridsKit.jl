# Multi-fidelity Modelling
Differentiation is supported for scalar valued functions currently.
This is achieved using automatic differentiation from the package `ReverseDiff.jl`.

```@example diff
using SparseGridsKit
n = 2
p(x) = [3*x[1]^3*x[2]^2]
p_prime(x1,x2) = [9*x1^2*x2^2, 6*x1^3*x2]

sg = create_sparsegrid(create_smolyak_miset(n, 4))
f_on_Z = p.(get_grid_points(sg))
f_sga = SparseGridApproximation(sg, f_on_Z)
f_sg_diff = derivative(f_sga)

# Create test points
test_points = [[x,y] for x in range(-1, stop=1, length=10), for y in range(-1, stop=1, length=10)]
# Evaluate the original function and its derivative at the test points
f_test = p_prime.(test_points)
f_test_diff = f_sg_diff.(test_points)
[f_test, f_test_diff, f_test - f_test_diff]
```

## Function Reference
```@autodocs
Modules = [SparseGridsKit]
Pages   = ["differentiation.jl"]
```