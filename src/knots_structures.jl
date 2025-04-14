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
    UniformPoints(domain::Vector{Real})
Generate uniform points on interval [domain[1],domain[2]].
# Arguments
- `domain` : Domain vector [a,b]
# Returns
- `UniformPoints` struct
"""
struct UniformPoints{T<:Real} <: Points
    domain::Vector{T}
end
UniformPoints() = UniformPoints([-1.0, 1.0])

"""
    (p::UniformPoints)(n)
Generate uniform points on domain.
# Arguments
- `n` : Number of points
"""
function (p::UniformPoints)(n)
    transformdomain(uniformpoints(n), p.domain[1], p.domain[2])
end

"""
    CCPoints(domain::Vector{Real})
Generate Clenshaw--Curtis points on interval [domain[1] domain[2]].
# Arguments
- `domain` : Domain
# Returns
- `CCPoints` struct
"""
struct CCPoints{T<:Real} <: Points
    domain::Vector{T}
end
CCPoints() = CCPoints([-1.0,1.0])

"""
    (p::CCPoints)(n)
Generate Clenshaw-Curtis points on domain.
# Arguments
- `n` : Number of points
"""
function (p::CCPoints)(n)
    transformdomain(ccpoints(n),p.domain[1], p.domain[2])
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
    domain::Vector{T}
end
GaussHermitePoints() = GaussHermitePoints([-Inf, Inf])


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
    GaussLegendrePoints(domain::Vector{Real})
Generate Gauss--Legendre points on interval [domain[1],domain[2]].
# Arguments
- `domain` : Domain
# Returns
- `GaussLegendrePoints` struct
"""
struct GaussLegendrePoints{T<:Real} <: Points
    domain::Vector{T}
end
GaussLegendrePoints() = GaussLegendrePoints([-1.0,1.0])
"""
    (p::GaussLegendrePoints)(n)
Generate Gauss-Legendre points on domain.
# Arguments
- `n` : Number of points
"""
function (p::GaussLegendrePoints)(n)
    transformdomain(gausslegendrepoints(n), p.domain[1], p.domain[2])
end

"""
    LejaPoints([a,b],type,v,symmetric)
    LejaPoints() = LejaPoints([-1,1], :optim, z->sqrt(0.5), false)
    LejaPoints(domain) = LejaPoints(domain, :optim, z->sqrt(0.5), false)
    LejaPoints(domain,type) = LejaPoints(domain,type, z->sqrt(0.5), false)


Generate Leja points on interval [domain[1],domain[2].
# Arguments
- `domain` : Domain
- `symmetric` : Boolean indicating if the points are symmetric (default is false).
- `type` : Type of Leja points (:optim or :classic), to use optimisation based or discrete search minimisation.
- `v` : Function to compute the weight function (default is sqrt(0.5)).
# Returns
- `LejaPoints` struct
"""
struct LejaPoints{T<:Real} <: Points
    domain::Vector{T}
    symmetric::Bool
    type::Symbol
    v::Function
end
LejaPoints() = LejaPoints([-1.0, 1.0], false, :optim, z->sqrt(0.5))
LejaPoints(domain) = LejaPoints(domain, false, :optim, z->sqrt(0.5))
LejaPoints(domain,symmetric) = LejaPoints(domain, symmetric, :optim, z->sqrt(0.5))

"""
    (p::LejaPoints)(n)
Generate Leja points on interval [a,b].
# Arguments
- `n` : Number of points
"""
function (p::LejaPoints)(n)
    if p.type == :optim
        xw = lejapoints(n; v=p.v, symmetric=p.symmetric)
    elseif p.type == :classic
        xw = lejapointsdiscretesearch(n, type=p.symmetric)
    else
        error("Unknown Leja points type: $(p.type)")
    end
    return transformdomain(xw, p.domain[1], p.domain[2])
end

"""
    CustomPoints(domain::Vector{Real}, knotsfunction::Function)
Generate custom points on interval [domain[1],domain[2]].

# Arguments
- `domain` : Domain
- `knotsfunction` : Function to generate knots from `n` points

# Returns
- `CustomPoints` struct
"""
struct CustomPoints{T<:Real} <: Points
    domain::Vector{T}
    knotsfunction::Function
end
"""
    (p::CustomPoints)(n)
Generate custom points on domain.

# Arguments
- `n` : Number of points
"""
function (p::CustomPoints)(n)
    p.knotsfunction(n)
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

struct CustomLevel <: Level
    f::Function
end
"""
    (level::CustomLevel)(l)
Custom rule
# Arguments
- `l` : Level number
- `function` : Function to generate knots from `n` points
# Returns
- ``n` : Number of knots
"""
function (level::CustomLevel)(l)
    return level.f(l)
end