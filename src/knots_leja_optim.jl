"""
    lejapoints(n; v=z->sqrt(0.5), method=GradientDescent(), symmetric=false)

Generate Leja quadrature points and weights on domain [a,b].

# Arguments
- `n`: Number of points.
- `v`: (Optional) Weight function. Usually `z -> sqrt(\rho(y))` where `\rho(y)`is the density function. Default `\rho=0.5`.
- `method`: (Optional) Optimisation method. Default is GradientDescent().
- `symmetric`: (Optional) If true, the points are symmetric about 0. Default is false.

# Returns
- `x`: Knots
- `w`: Weights
"""
function lejapoints(n; v=z->sqrt(0.5), method=GradientDescent(), symmetric=false)    
    leja = Vector{Float64}(undef,n)
    w = Vector{Float64}(undef,n)
    leja[1] = 1.0

    # Use Chebyshev candidates for approximate minimiser
    M = 1e2+1 # Number of candidates
    X = cos.((0:(M-1)) .* π / (M-1))

    f(z, Z) = - v(z)^2 * prod((z .- Z).^2)

    # Use symm_count to create symmetry if necessary, whilst ignoring repeating 0.0.
    symm_count = 2
    for kk in 2:n
        # Enforce symmetry if required
        if iseven(symm_count) && symmetric
            leja[kk] = -leja[kk-1]
        else
            g(z) = f(z,leja[1:kk-1])

            # Using the constrained search did not give the expected results.
            # The sixth point was different to previous results.
            # results = optimize(g,-1.0,1.0, Brent())
            # For this reason, a coarser search via candidates is used to find an initial guess.
            # GradientDescent() is then applied.
            # Approximate minimizer
            initial_x = X[argmax(v.(X).*prod(abs.(X .- leja[1:(kk-1)]'), dims=2))]
            results = optimize(g, [initial_x], method, Optim.Options(f_abstol=1e-14); autodiff= :forward)
            leja[kk] = results.minimizer[1]
        end
        if !isapprox(leja[kk],0.0;atol=1e-14)
            symm_count += 1
        end
    end

    # Compute weights
    p = Vector{Fun}(undef,1)
    p_eval = Vector{Float64}(undef, n)
    (x_cc,w_cc) = CCPoints()(n)
    for ii in eachindex(leja)
        p[1] = poly_Fun(leja, ii)
        lagrange_evaluation!(p_eval,leja,ii,x_cc)
        w[ii] = (v(x_cc).^2 .* w_cc)' * p_eval
    end
    
    return leja, w
end