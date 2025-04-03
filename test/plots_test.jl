using Plots
@testset "Sparse Grid Adaptivity Tests" begin

    p = plot(CCPoints(-3,3))
    try
        plot!(CCPoints(-3,3); npts=21)
        plotpoints == true
    catch
        plotpoints = false
    end
    @test plotpoints

    miset = create_smolyak_miset(3,3)
    try
        p = misetplot(miset)
        plotmiset = true
    catch
        plotmiset = false
    end
    @test plotmiset


    sg = create_sparsegrid(miset,fill([-1,1],3))
    p = plot(sg; targetdims=[3,2,1])

    f(x) = x[1]^5 + cos(x[2]) + abs(x[3])
    f_on_sg = f.(get_grid_points(sg))

    sga = SparseGridApproximation(sg,f_on_sg)
    ssg = convert_to_spectral_approximation(sga)

    try
        p = plot(
            plot(sga; fill=true),
            plot(sga; targetdims=[2,3]),
            plot(ssg; color=:turbo, fill=true),
            plot(ssg; seriestype=:surface, targetdims=[2,3]),
            layout = (2,2)
        )
        plotapproximations = true
    catch
        plotapproximations = false
    end
    @test plotapproximations

end