abstract type Points end

"""
    Points(n::Integer)
Generate knots on interval [a,b] using the specified generator function.
# Arguments
- `n` : Number of points
"""
function (p::Points)(n)
    error("$(typeof(p)) does not implement `(::$(typeof(p)))(n)`")
end

"""
    UniformPoints(a::Real, b::Real)
Generate uniform points on interval [a,b].
# Arguments
- `a` : Lower bound
- `b` : Upper bound
# Returns
- `UniformPoints` struct
"""
struct UniformPoints{T<:Real} <: Points
    a::T
    b::T
end
UniformPoints() = UniformPoints(-1.0, 1.0)

"""
    (p::UniformPoints)(n)
Generate uniform points on interval [a,b].
# Arguments
- `n` : Number of points
"""
function (p::UniformPoints)(n)
    uniformpoints(n, p.a, p.b)
end

"""
    CCPoints(a::Real, b::Real)
Generate Clenshaw--Curtis points on interval [a,b].
# Arguments
- `a` : Lower bound
- `b` : Upper bound
# Returns
- `CCPoints` struct
"""
struct CCPoints{T<:Real} <: Points
    a::T
    b::T
end
CCPoints() = CCPoints(-1.0,1.0)

"""
    (p::CCPoints)(n)
Generate Clenshaw-Curtis points on interval [a,b].
# Arguments
- `n` : Number of points
"""
function (p::CCPoints)(n)
    ccpoints(n, p.a, p.b)
end

"""
    GaussHermitePoints()
Generate Gauss--Hermite points.
# Arguments
- `None`
# Returns
- `GaussHermitePoints` struct
"""
struct GaussHermitePoints{T<:Real} <: Points
end
GaussHermitePoints() = GaussHermitePoints{Float64}()


"""
    (p::GaussHermitePoints)(n)
Generate Gauss-Hermite points on interval (-∞,∞).
# Arguments
- `n` : Number of points
"""
function (p::GaussHermitePoints)(n)
    gausshermitepoints(n)
end

"""
    GaussLegendrePoints(a::Real, b::Real)
Generate Gauss--Legendre points on interval [a,b].
# Arguments
- `a` : Lower bound
- `b` : Upper bound
# Returns
- `GaussLegendrePoints` struct
"""
struct GaussLegendrePoints{T<:Real} <: Points
    a::T
    b::T
end
GaussLegendrePoints() = GaussLegendrePoints(-1.0,1.0)
"""
    (p::GaussLegendrePoints)(n)
Generate Gauss-Legendre points on interval [a,b].
# Arguments
- `n` : Number of points
"""
function (p::GaussLegendrePoints)(n)
    gausslegendrepoints(n, p.a, p.b)
end

"""
    LejaPoints(a::Real, b::Real)
Generate Leja points on interval [a,b].
# Arguments
- `a` : Lower bound
- `b` : Upper bound
# Returns
- `LejaPoints` struct
"""

struct LejaPoints{T<:Real} <: Points
    a::T
    b::T
end
LejaPoints() = LejaPoints(-1.0,1.0)
"""
    (p::LejaPoints)(n)
Generate Leja points on interval [a,b].
# Arguments
- `n` : Number of points
"""
function (p::LejaPoints)(n)
    lejapoints(n, p.a, p.b)
end

abstract type Level end

"""
    (level::Level)(n::Integer)
Level to knot number mapping
# Arguments
- `l` : Level number
# Returns
- ``n` : Number of knots
"""
function (level::Level)(l)
    error("$(typeof(level)) does not implement `(::$(typeof(level)))(n)`")
end

struct Doubling <: Level
end
"""
    (level::Doubling)(l)
Doubling rule
# Arguments
- `l` : Level number
# Returns
- ``n` : Number of knots
"""
function (level::Doubling)(l)
    return doubling(l)
end

struct Tripling <: Level
end
"""
    (level::Tripling)(l)
Tripling rule
# Arguments
- `l` : Level number
# Returns
- ``n` : Number of knots
"""
function (level::Tripling)(l)
    return tripling(l)
end

struct Linear <: Level
end
"""
    (level::Linear)(l)
Linear rule
# Arguments
- `l` : Level number
# Returns
- ``n` : Number of knots
"""

function (level::Linear)(l)
    return linear(l)
end

struct TwoStep <: Level
end
"""
    (level::TwoStep)(l)
TwoStep rule
# Arguments
- `l` : Level number
# Returns
- ``n` : Number of knots
"""
function (level::TwoStep)(l)
    return twostep(l)
end