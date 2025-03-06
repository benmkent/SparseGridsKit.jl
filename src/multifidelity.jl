"""
    fidelity(n)

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
function multifidelityfunctionwrapper(f, nfid, nparam)
        return z -> f(z[1:nfid],z[nfid+1:-1])
end