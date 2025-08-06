using Plots
@testset "Sparse Grid Adaptivity Tests" begin
    # This test set is simple, it will fail if the plotting fails but otherwise tests nothing
    p = plot(CCPoints([-3,3]))
    plot!(CCPoints([-3,3]); npts=21)
    
    @test true

    sg = create_sparsegrid(miset)
    p = plot(sg; targetdims=[3,2,1])

    f(x) = x[1]^5 + cos(x[2]) + abs(x[3])
    f_on_sg = f.(get_grid_points(sg))

    sga = SparseGridApproximation(sg,f_on_sg)
    ssg = convert_to_spectral_approximation(sga)

    p = plot(
        plot(sga; fill=true),
        plot(sga; targetdims=[2,3]),
        plot(ssg; color=:turbo, fill=true),
        plot(ssg; seriestype=:surface, targetdims=[2,3]),
        layout = (2,2)
    )
    @test true

end