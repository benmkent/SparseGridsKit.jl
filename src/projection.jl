using Optim

function project_onto_sparsegrid(sg::SparseGrid, x::Vector, y::Vector)

    # Construct Vandermode matrix
    n = length(x)
    d = length(x[1])
    n_grid_pts = get_n_grid_points(sg)
    A = zeros(n, n_grid_pts)
    for i in 1:n_grid_pts
        f_on_sg = zeros(n_grid_pts)
        f_on_sg[i] = one(1) 
        A[:,i] = interpolate_on_sparsegrid(sg, f_on_sg, x)
    end

    # Solve least squares problem
    f_optim = (u) -> sum((y - A * u).^2) + sum(convert_to_spectral_approximation(SparseGridApproximation(sg,u)).coefficients.^2)

    
    sga = SparseGridApproximation(sg, f_lsq_on_Z)
    return sga
end