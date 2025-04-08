@testset "Custom knots" begin
    cp = CustomPoints([1.0,10.0], n->(1:n, 1/n*ones(n)))
    cl = CustomLevel(l->l)
    x,w = cp(5)
    @test x ≈ 1:5
    @test sum(w) ≈ 1.0
end

@testset "Leja points" begin
    # Test Leja points
    n = 5
    lp = LejaPoints()
    x,w = lp(n)
    @test length(x) == n
    @test length(w) == n
    @test minimum(x) ≈ -1.0 && maximum(x) ≈ 1.0
    @test sum(w) ≈ 1.0

    symmetric = true
    lp = LejaPoints([-1.0,1.0],symmetric) # symmetric points
    x,w = lp(21)
    sort!(x)
    @test all(x[ii] ≈ -x[end-ii+1] for ii in 1:11)
    @test sum(w) ≈ 1.0
end

@testset "Mixed Grids Tests" begin
    # Test create_sparsegrid
    n,k =3,3
    knots = [GaussHermitePoints(), CCPoints([-2.0,10.0]), UniformPoints([0.0,1.0])]
    rules = [Linear(), Doubling(), Doubling()]
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set,knots=knots, rule=rules)

    Z = hcat(get_grid_points(sg)...)'

    Z_cc = Z[:,2]
    Z_uniform = Z[:,3]

    @test minimum(Z_cc) ≈ -2.0
    @test maximum(Z_cc) ≈ 10.0
    @test minimum(Z_uniform) ≈ 0.0
    @test maximum(Z_uniform) ≈ 1.0
end