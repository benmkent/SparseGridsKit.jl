# Sparse Grids for More General Objects
More general objects can be interpolated and integrated via `SparseGridsKit.jl`.
The SparseGridsKit supports objects for which addition is defined.
For example, we can could consider [ApproxFun](https://github.com/JuliaApproximation/ApproxFun.jl) objects or even Python objects provided they support addition.

## Example: ApproxFun.jl
### One Dimensional Functions
We can consider interpolation of `Fun` objects.
First a parametric `Fun` object can be evaluated at sparse grid points.
```@example approxfun
using SparseGridsKit, ApproxFun, Plots, LaTeXStrings
n, k = 1, 1
mi_set = create_smolyak_miset(n, k)
sg = create_sparsegrid(mi_set)
f = y-> Fun(Taylor(),[y,2,3])
f_on_grid = [f(x[1]) for x in get_grid_points(sg)]
```
The `Fun` objects can then be interpolated using the `interpolate_on_sparsegrid` function.
The interpolated `Fun` objects are plotted below.
```@example approxfun
target_points = [[x] for x in -1.0:0.2:1.0]
interpolation_result = interpolate_on_sparsegrid(sg, f_on_grid, target_points)
plot(interpolation_result,
    xlims = (-2,2),
)
```
Integration is also possible using the sparse grid formulation.
```@example approxfun
pcl = precompute_lagrange_integrals(4)
expectedvalue = integrate_on_sparsegrid(sg,f_on_grid,pcl)
#plot(expectedvalue)
```
### Two dimensional Elliptic PDE
The `ApproxFun` package offers a fast way to solve a model elliptic PDE.
Given that the `Fun` objects can be used naturally within the `SparseGridsKit.jl`, we can consider computing an approximation to the solution of a parameteric elliptic PDE.

This is set up on a two dimensional spatial domain $D=[-1,1]^2$.
```@example approxfun
d = ChebyshevInterval()^2
```
A parametric forcing function is defined with a parameter domain $\Gamma=[-1,1]^2$.
```@example approxfun
x,y = Fun(d)
# Forcing is parametric
f = z -> ProductFun((x,y)->exp.(-(1.01+z[1])*(x+0.3)^2-2(y-z[2])^2))  # use broadcasting as exp(f) not implemented in 2D
```
The forced elliptic PDE
```math
-\nabla u(x,z) = f(x,z) \text{ for all } x \in D, u(x,z) = 0 \text{ for all } x \in \partial D
```
is easily solved pointwise for any parameter $z \in \Gamma$.
```@example approxfun
A = [Dirichlet(d); -Laplacian()]
# Sln is parametric
u = z -> ProductFun(\(A, [zeros(∂(d)); f(z)]; tolerance=1E-6),d)

z = [0, 0]
uz = u(z)
plot(uz)
plot!(xlabel=L"x_1",
        ylabel=L"x_2",
        zlabel=L"u(x,z)")
```

A sparse grid approximation can be constructed.
The pointwise solutions on the sparse grid can be interpolated to give an approximation at any point in the domain.
The animation below illustrates the changing sparse grid approximation and the approximation error for parameters on the circle $\Vert x \Vert_2 = 1$.
```@example approxfun
n, k = 2,2
mi_set = create_smolyak_miset(n, k)
sg = create_sparsegrid(mi_set)
f_on_grid = [u(z) for z in get_grid_points(sg)]

nsteps = 20
points = get_grid_points(sg)
x = [p[1] for p in points]
y = [p[2] for p in points]
@gif for i in range(0, stop = 2π, length = nsteps)
        z = [cos(i),sin(i)]
        interpolation_result = interpolate_on_sparsegrid(sg, f_on_grid, [z])

        p1 = plot(interpolation_result[1])
        plot!(zlims!(0,0.3),title="SG Approximation")
        plot!(xlabel=L"x_1",
        ylabel=L"x_2",
        zlabel=L"u(x,z)")

        p2 = plot(interpolation_result[1]-u(z))
        plot!(zlims!(-3e-2,3e-2),title="Interpolation Error")
        plot!(xlabel=L"x_1",
        ylabel=L"x_2",
        zlabel=L"e(x,z)")

        p3 = scatter(x,y,label="grid")
        scatter!(p3, [z[1]],[z[2]],label="z")
        plot(p1, p2, p3, layout=[2,2])
        plot!(xlabel=L"y_1",
        ylabel=L"y_2")
end
```
The sparse grid has been constructed using Clenshaw--Curtis points.
The quadrature rules can be applied assuming an underlying weight function $\rho(y)=0.5^2$.
If this is interpreted as a probabiliy then integration gives the expected value of the PDE solution.
```@example approxfun
pcl = precompute_lagrange_integrals(k+1)

expectedvalue = ProductFun(integrate_on_sparsegrid(sg,f_on_grid, pcl),d)
plot(expectedvalue)
plot!(title="Approximation of Expected Value",
        xlabel=L"x_1",
        ylabel=L"x_2",
        zlabel=L"\mathbb{E}[u(x,\cdot)]")
```
Similarly an approximation of the variance can be computed and interpolated.
```@example approxfun
square = PF -> ProductFun(PF*PF,d)
variance = ProductFun(integrate_on_sparsegrid(sg,[square(f_on_grid[i]-expectedvalue) for i in eachindex(f_on_grid)], pcl),d)
plot(variance)
plot!(title="Approximation of Variance",
        xlabel=L"x_1",
        ylabel=L"x_2",
        zlabel=L"Var[u(x,\cdot)]")
```