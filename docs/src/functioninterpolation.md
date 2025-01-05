# Interpolation of PDE Solution Objects
More general objects can be interpolated.
The SparseGridsKit supports objects for which addition is defined.
For example, we can could consider [ApproxFun](https://github.com/JuliaApproximation/ApproxFun.jl) objects.
Or even Python objects.

## ApproxFun.jl Funs
### One Dimensional Functions
We can consider interpolation of `Fun` objects.
```@example approxfun
using SparseGridsKit, ApproxFun, Plots
n, k = 1, 1
mi_set = create_smolyak_miset(n, k)
sg = create_sparsegrid(mi_set)
f = y-> Fun(Taylor(),[y,2,3])
f_on_grid = [f(x[1]) for x in get_grid_points(sg)]
```
The `Fun` objects can then be interpolated.
```@example approxfun
target_points = [[x] for x in -1.0:0.2:1.0]
interpolation_result = interpolate_on_sparsegrid(sg, f_on_grid, target_points)
plot(interpolation_result)
```
Integration is also possible.
This can be interpreted as computing the expected value.
```@example approxfun
pcl = precompute_lagrange_integrals(4)
expectedvalue = integrate_on_sparsegrid(sg,f_on_grid,pcl)
#plot(expectedvalue)
```
### Two dimensional Elliptic PDE
The `ApproxFun` package offers a fast way to solve a model elliptic PDE.
This is set up on a two dimensional spatial domain $D=[-1,1]^2$.
```@example approxfun
d = ChebyshevInterval()^2
```
The forcing is defined to be parametric taking two variables.
```@example approxfun
x,y = Fun(d)
# Forcing is parametric
f = z -> ProductFun((x,y)->exp.(-(1.01+z[1])*(x+0.3)^2-2(y-z[2])^2))  # use broadcasting as exp(f) not implemented in 2D
```
This is then easily solved for any parameter realisation.
```@example approxfun
A = [Dirichlet(d); -Laplacian()]
# Sln is parametric
u = z -> ProductFun(\(A, [zeros(∂(d)); f(z)]; tolerance=1E-6),d)

z = [0, 0]
uz = u(z)
plot(uz)
```

A sparse grid approximation is then constructed.
This can be interpolated to give an approximation at 
```@example approxfun
n, k = 2,2
mi_set = create_smolyak_miset(n, k)
sg = create_sparsegrid(mi_set)
f_on_grid = [u(z) for z in get_grid_points(sg)]

nsteps = 100
points = get_grid_points(sg)
x = [p[1] for p in points]
y = [p[2] for p in points]
@gif for i in range(0, stop = 2π, length = nsteps)
        z = [cos(i),sin(i)]
        interpolation_result = interpolate_on_sparsegrid(sg, f_on_grid, [z])
        p1 = plot(interpolation_result[1])
        plot!(zlims!(0,0.3),title="SG Approximation")
        p2 = plot(interpolation_result[1]-u(z))
        plot!(zlims!(-3e-2,3e-2),title="Interpolation Error")
        p3 = scatter(x,y,label="grid")
        scatter!(p3, [z[1]],[z[2]],label="z")
        plot(p1, p2, p3, layout=[2,2])
end
```
Integration can then be used to give the expected value and the variance.
```@example approxfun
pcl = precompute_lagrange_integrals(k+1)

expectedvalue = ProductFun(integrate_on_sparsegrid(sg,f_on_grid, pcl),d)
plot(expectedvalue)
```
```@example approxfun
square = PF -> ProductFun(PF*PF,d)
variance = ProductFun(integrate_on_sparsegrid(sg,[square(f_on_grid[i]-expectedvalue) for i in eachindex(f_on_grid)], pcl),d)
plot(variance)
```