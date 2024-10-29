"""
Generates the Genz test functions.

# Arguments
- `n::Int`: Number of dimensions (positive integer).
- `C::Float64`: Scaling coefficient.
- `W::Float64`: Shift parameter.
- `T::String`: Decay type, options are "nodecay", "quadraticdecay", "quarticdecay", "exponentialdecay", and "squaredexponentialdecay".
- `N::String`: Function type, options are "oscillatory", "productpeak", "productpeaknormalised", "cornerpeak", "gaussianpeak", "c0continuous", and "discontinuous".

# Returns
- A function of specified Genz test function.
"""
function genz(n::Int, C::Float64, W::Float64, T::String, N::String)
    # Minimum coefficient
    cmin = 5e-6

    # Generate coefficients based on the decay type
    c = Vector{Float64}(undef, n)
    if T == "nodecay"
        c .= (1:n .+ 0.5) / n
    elseif T == "quadraticdecay"
        c .= 1.0 ./ ((1:n) .+ 1).^2
    elseif T == "quarticdecay"
        c .= 1.0 ./ ((1:n) .+ 1).^4
    elseif T == "exponentialdecay"
        c .= exp.(log(cmin) * ((1:n) .+ 1) ./ n)
    elseif T == "squaredexponentialdecay"
        c .= 10.0 .^ (log10(cmin) * ((1:n) .+ 1) / n)
    else
        error("unknown coefficients")
    end

    # Normalize coefficients and apply scaling
    c = C * c / sum(c)

    # Initialize weights
    w = fill(W, n)

    # Define the function based on the specified type
    if N == "oscillatory"
        return θ -> cos(2π * w[1] + sum(c .* θ))
    elseif N == "productpeak"
        return θ -> prod(c.^(-2) .+ (θ .- w).^2)^(-1)
    elseif N == "productpeaknormalised"
        fpp = genz(n, C, W, T, "productpeak")
        return θ -> fpp(θ) * prod(c.^(-2))
    elseif N == "cornerpeak"
        return θ -> (1 + sum(c .* θ))^(-(n + 1))
    elseif N == "gaussianpeak"
        return θ -> exp(-sum(c.^2 .* (θ .- w).^2))
    elseif N == "c0continuous"
        return θ -> exp(-sum(c .* abs.(θ .- w)))
    elseif N == "discontinuous"
        return θ -> (θ[1] <= w[1]) * (θ[2] <= w[2]) * exp(sum(c .* θ))
    else
        error("unknown integrand family")
    end
end
