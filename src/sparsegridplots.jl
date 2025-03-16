export plot_pairwise_mi, plot_sparse_grid

using Plots, StatsPlots
"""
    plot_pairwise_mi(miset)
Plots the pairwise histograms for each variable pair.
"""
function plot_pairwise_mi(miset)
    mi_vec = get_mi(miset)
    mi_vec_matrix=transpose(hcat(mi_vec...))
    nmi, ndim = size(mi_vec_matrix)
    psub = Matrix{Any}(undef,ndim,ndim)
    for ii=1:ndim, jj=1:ndim
        countdict =StatsPlots.countmap([mi_vec_matrix[kk,[ii,jj]] for kk=1:nmi])
        keyvalues = collect(keys(countdict))
        cvalues = collect(values(countdict))
        x = [k[1] for k in keyvalues]
        y = [k[2] for k in keyvalues]
        z = [c[1] for c in cvalues]
        psub[ii,jj] = scatter(x,y, zcolor=collect(values(countdict)),legend=:none)
        if ii == 1
        plot!(ylabel="Parameter "*string(jj))
        end
        if jj == 1
        plot!(xlabel="Parameter "*string(ii))
        end
    end
    p=plot([psub[ii,jj] for jj=ndim:-1:1 for ii=1:ndim]...,layout=(ndim,ndim))
end

"""
    plot_sparse_grid(sg)
Plots the sparse grid points
# Arguments
- `sg`: Sparse grid to plot
- `dims`: (optional) for sparse grids with dimension greater than 3, then selects which 3 dimensions to plot as xyz.

"""
function plot_sparse_grid(sg; dims=nothing)
    pts = get_grid_points(sg)
    p = plot()

    if sg.dims == 1
        scatter!(p, [pt[1] for pt in pts], fill(1, length(pts)))
    elseif sg.dims == 2
        scatter!(p, [pt[1] for pt in pts], [pt[2] for pt in pts])
    elseif isnothing(dims)
        scatter!(p, [pt[1] for pt in pts], [pt[2] for pt in pts], [pt[3] for pt in pts])
    elseif length(dims) == 3
        scatter!(p, [pt[dims[1]] for pt in pts], [pt[dims[2]] for pt in pts], [pt[dims[3]] for pt in pts])
    else
        error("Invalid dims argument. Provide exactly 3 dimensions for higher-dimensional plots.")
    end

    return p
end

