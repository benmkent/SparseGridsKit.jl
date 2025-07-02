import PolyChaos: clenshaw_curtis
import FastGaussQuadrature: gausslegendre, gausshermite



"""
    transformdomain(xw, a, b)

# Arguments
- `xw`: Points and weights tuple
- `a`, `b`: Domain [a,b]

# Returns
- Transformed x and w.
"""
transformdomain(xw, a, b) = (a .+ (b-a) .* 0.5 .* (xw[1] .+ 1.0), xw[2])

"""
    doubling(n)

Use doubling rule for levels `2^{n-1} + 1` for `n > 1`.

# Arguments
- `n`: Level

# Returns
- Number of points
"""
doubling(n) = n == 1 ? 1 : 2^(n - 1) + 1

"""
    tripling(n)

Use tripling rule for levels `3^{n-1} + 1` for `n > 1`.

# Arguments
- `n`: Level

# Returns
- Number of points
"""
tripling(n) = n == 1 ? 1 : 3^(n - 1)

"""
    linear(n)

Use linear level to points

# Arguments
- `n`: Level.

# Returns
- Number of points
"""
linear(n) = n

"""
    twostep(n)

Use twostep level to points

# Arguments
- `n`: Level.

# Returns
- Number of points
"""
twostep(n) = 2*(n-1)+1

"""
    uniformpoints(n)

Uniformly spaced quadrature points on interval [-1,1]

# Arguments
- `n`: Number of points.

# Returns
- Points and weights in `[-1, 1]`.
"""
uniformpoints(n) = n == 1 ? ([0.0], [1.0]) : (range(-1, stop=1, length=n), 1/n .* ones(n))


"""
    ccpoints(n)

Generate Clenshaw-Curtis quadrature points and weights on domain `[-1, 1]`.

# Arguments
- `n`: Number of points.

# Returns
- Points and weights in `[-1, 1]`.
"""
function ccpoints(n)
    if n == 1
        x = [0.0]
        w = [1.0]
    elseif n == 2
        x = [-1.0, 1.0]
        w = [0.5, 0.5]
    else
        x, w = clenshaw_curtis(n - 1)
        x = 0.5 * (x .- x[end:-1:1]) # Force symmetry
        w = 0.5 * w
    end
    return x, w
end

"""
    gausslegendrepoints(n)

Generate Gauss-Legendre quadrature points and weights on domain `[-1, 1]`.

# Arguments
- `n`: Number of points.

# Returns
- Points and weights in `[-1, 1]`.
"""
function gausslegendrepoints(n)
    x,w = gausslegendre(n);
    x = x
    w = w/2.0
    return x, w
end


"""
    gausshermitepoints(n)

Generate Gauss-Hermite quadrature points and weights on domain `[-∞, ∞]`.

# Arguments
- `n`: Number of points.

# Returns
- Points and weights in `[-∞, ∞]`.
"""
function gausshermitepoints(n)
    x,w = gausshermite(n);
    x=x*sqrt(2)
    w = w/sqrt(π)
    return x, w
end

"""
    lejapointsdiscretesearch(n)

Generate Leja quadrature points and weights on domain `[-1, 1]` using a search over Chebyshev candidates.

# Arguments
- `n`: Number of points.

# Returns
- `x`: Knots
- `w`: Weights
"""
function lejapointsdiscretesearch(n; type=nothing)
    # Use Chebyshev candidates
    M = 1e3+1 # Number of candidates
    X = cos.((0:(M-1)) .* π / (M-1))
    
    leja = Vector{Float64}(undef,n)
    w = Vector{Float64}(undef,n)
    leja[1] = X[argmax(abs.(X))]

    for kk in 2:n
        if kk > 3 && isodd(kk) && type == :symmetric
            next_x = -leja[kk-1]
        else
            next_x = X[argmax(prod(abs.(X .- leja[1:(kk-1)]'), dims=2))]
        end
        leja[kk] = next_x
    end

    # Compute weights
    p = Vector{Fun}(undef,1)
    for ii in eachindex(leja)
        p[1] = poly_Fun(leja, ii)
        w[ii] = sum(p[1])*0.5
    end
    
    return leja, w
end

"""
    lejapoints(n; v=z->sqrt(0.5), method=GradientDescent(), symmetric=false)

Generate Leja quadrature points and weights on domain [-1,1].

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
    for ii in eachindex(leja)
        p[1] = poly_Fun(leja, ii)
        w[ii] = sum(p[1])*0.5
    end
    
    return leja, w
end
