using DifferentiationInterface
backend = AutoForwardDiff()

"""
    derivative(sga::SparseGridApproximation)
Computes the derivative of a SparseGridApproximation object.
# Arguments
- `sga`: An instance of `SparseGridApproximation`.
# Returns
-`SparseGridApproximation` object representing the derivative.
"""
function derivative(sga::SparseGridApproximation)
    ssg = convert_to_spectral_approximation(sga)
    ssg_diff = derivative(ssg; sparsegrid=sga.sg)
    # Convert back to sparse grid approximation
    grid = get_grid_points(sga.sg)
    f_diff_on_grid = ssg_diff.(grid)
    # Create a new SparseGridApproximation object
    sga_diff = SparseGridApproximation(sga.sg, f_diff_on_grid)
    return sga_diff
end
"""
    derivative(sg::SparseGrid, f_on_z::Vector)
Computes the derivative of a function defined on a sparse grid.
# Arguments
- `sg`: A sparse grid object.
- `f_on_z`: A vector of function values on the sparse grid.
# Returns
- A tuple containing the sparse grid and the derivative values on the grid.
"""
function derivative(sg::SparseGrid, f_on_z::Vector)
    sga = SparseGridApproximation(sg,f_on_z)
    ssg_diff = derivative(sga)
    return ssg_diff.sg, ssg_diff.fongrid
end
"""
    derivative(ssg::SpectralSparseGridApproximation; sparsegrid=nothing)
Computes the derivative of a SpectralSparseGridApproximation object.
# Arguments
- `ssg`: An instance of `SpectralSparseGridApproximation`.
- `sparsegrid`: Optional sparse grid to use for differentiation.
# Returns
- A `SpectralSparseGridApproximation` object representing the derivative.
"""
function derivative(ssg::SpectralSparseGridApproximation; sparsegrid=nothing)
    # Get corresponding sparse grid
    if isnothing(sparsegrid)
        knots = Vector{Points}(undef, length(ssg.polytypes))
        rules = Vector{Level}(undef, length(ssg.polytypes))
        # Decide on knots and rules
        for (ii,pt) in enumerate(ssg.polytypes)
            if isa(pt,Chebyshev)
                knots[ii] = CCPoints(getdomain(pt))
                rules[ii] = Doubling()
            elseif isa(pt,Hermite)
                knots[ii] = GaussHermitePoints()
                rules[ii] = Doubling()
            else
                @error "Unsupported underlying polynomial type" pt
                error()
            end
        end
        # Convert spectral approximation basis to sparse grid
        sga = convert_to_sg_approximation(ssg, knots, rules)
        sparsegrid = sga.sg
    end
    Z = get_grid_points(sparsegrid)
    grad_on_Z = Vector(undef, length(Z))
    for (ii,z) in enumerate(Z)

        grad_on_Z[ii] = jacobian(y->ssg(y), backend, z)
    end
    sga_diff = SparseGridApproximation(sparsegrid,grad_on_Z)
    ssg_diff = convert_to_spectral_approximation(sga_diff)
    return ssg_diff
end