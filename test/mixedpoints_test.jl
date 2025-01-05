@testset "Mixed Grids Tests" begin
    # Test create_sparsegrid
    n,k =3,0
    knots = [linearpoints, ccpoints, uniformpoints]
    rules = [linear, doubling, doubling]
    mi_set = create_smolyak_miset(n,k)
    sg = create_sparsegrid(mi_set,knots=knots, rule=rules)
end