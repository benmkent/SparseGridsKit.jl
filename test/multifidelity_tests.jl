using Test
using SparseGridsKit

@testset "Multi-fidelity Modelling" begin
    # FidelityPoints()(3), Fidelity(3)

    f(y) = 100*y[1] + 10*y[2] + 1*y[3]

    nfid = 1
    nparams = 3
    f_fidelities(alpha,y) = f(y) + 10.0^(-alpha[1])

    maxfidelity = 5
    rule = [fill(Fidelity(),nfid)..., fill(Doubling(),nparams)...]
    knots = [fill(FidelityPoints(),nfid)..., fill(CCPoints(),nparams)...]

    f_wrapped = multifidelityfunctionwrapper(f_fidelities,knots)

    MI = create_smolyak_miset(nfid+nparams,2)

    sg = create_sparsegrid(MI; knots=knots, rule=rule)
    f_eval = f_wrapped.(get_grid_points(sg))
    f_true = [10.0^-z[1] + 100.0*z[2]+10.0*z[3]+z[4] for z in get_grid_points(sg)]
    @test f_eval ≈ f_true

    pcl = precompute_lagrange_integrals(5, knots, rule)

    integral_result = integrate_on_sparsegrid(sg, f_eval, pcl)
    @test integral_result ≈ 0.001

    (sg, f_on_Z) = adaptive_sparsegrid(f_wrapped, nfid+nparams; maxpts = 100, proftol=1e-4, rule = rule, knots = knots, θ=1e-4, type=:deltaint, costfunction=α->10.0^prod(α))

    integral_result = integrate_on_sparsegrid(sg, f_on_Z, pcl)
    @test abs(integral_result) ≤ 1e-4 
end
