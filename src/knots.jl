import PolyChaos: clenshaw_curtis

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
    linear(n)

Use linear level to points

# Arguments
- `n`: Level.

# Returns
- Number of points
"""
linear(n) = n

"""
    linearpoints(n)

Generate points as Natural numbers 1,2,3,...

# Arguments
- `n`: Number of points.

# Returns
- `z` : 1:n
- `w` : ones(n)
"""
linearpoints(n) = 1:n, ones(n)

"""
    uniformpoints(n, a, b)

Uniformly spaced quadrature points on interval [a,b]

# Arguments
- `n`: Number of points.
- `a`, `b`: Interval bounds.

# Returns
- Points and weights in `[a, b]`.
"""
uniformpoints(n, a, b) = transformdomain(uniformpoints(n), a, b)

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
    ccpoints(n, a, b)

Generate Clenshaw-Curtis quadrature points and weights on domain [a,b].

# Arguments
- `n`: Number of points.
- `a`, `b`: Interval bounds.

# Returns
- Points and weights in `[a, b]`.
"""
ccpoints(n, a, b) = transformdomain(ccpoints(n), a, b)

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
    else
        x, w = clenshaw_curtis(n - 1)
        x = 0.5 * (x .- x[end:-1:1]) # Force symmetry
        w = 0.5 * w
    end
    return x, w
end
