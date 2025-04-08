struct FidelityPoints <: Points
    domain::Vector{Integer}
end
"""
    FidelityPoints()
Generate a struct for fidelity points
# Returns
- `FidelityPoints` struct
"""
FidelityPoints() = FidelityPoints([1, typemax(Int)])

"""
    (p::FidelityPoints)(n)
Generate fidelity points
# Arguments
- `n` : Number of points
"""
function (p::FidelityPoints)(n)
    fidelitypoints(n)
end

struct Fidelity <: Level
end
"""
    (f:Fidelity)(l)
Fidelity level to knot number mapping
# Arguments
- `l` : Fidelity level

# Returns
- `1`` : Fidelity level always returns 1`
"""
(f::Fidelity)(l) = fidelity(l)

"""
    fidelitypoints(n)

Generate indexes as Natural numbers 1,2,3,... for use as fidelity levels

# Arguments
- `n`: Fidelity index.

# Returns
- `z` : n
- `w` : 1.0
"""
fidelitypoints(n) = [n], [1.0]

fidelity(n) = 1

"""
    fidelitymap(y, nfid, nparam)

Splits a grid point y to y[1:nfid],y[nfid+1:nfid+nparam]

# Arguments
- `y`: Grid point
- `nfid`: Number of fidelity indices
- `nparam`: Number of parameters

# Returns
- `yfid` : fidelity indices
- `yparam` : parameters
"""
fidelitymap(y, nfid, nparam) = y[1:nfid], y[nfid+1:nfid+nparam]

"""
    multifidelityfunctionwrapper(f, nfid, nparam)

Wrap a function for multifidelity evaluation

# Arguments
- `f`: Function to wrap
- `z`: Grid point
- `nfid`: Number of fidelity indices
- `nparam`: Number of parameters
"""
function multifidelityfunctionwrapper(f, knots)
    if !isa(knots, Vector)
        f_wrapped = z -> f(z)
    else
        fidelity_mask = isa.(knots,FidelityPoints)
        f_wrapped = z -> f(z[fidelity_mask],z[.!fidelity_mask])
    end
    return f_wrapped
end