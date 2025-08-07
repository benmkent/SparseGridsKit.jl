using RecipesBase

"""
    f(p::Points; n = 11)

Recipe to plot scatter of knots and weights.

# Arguments
- `p::Points`: `Points` object.
- `n::Int`: (Optional) Number of knots. Defaults to `11`.
"""
@recipe function f(p::Points; n = 11)
    seriestype  :=  :scatter
    label --> string(typeof(p))*" "*string(n)
    xlabel --> "Knots"
    ylabel --> "Weight"
    x,w = p(n)
end

"""
    f(sg::SparseGrid; targetdims=[1,2,3])

    Recipe to convert sparse grid to set of x,y,z points in dimensions given by targetdims

# Arguments
- `sg` : Sparse grid
- `targetdims` : (Optional) Cross section dimensions
"""
@recipe function f(sg::SparseGrid; targetdims=[1,2,3])
    points = get_grid_points(sg)
    points_matrix = hcat(points...)'

    title --> "Sparse Grid"
    xlabel --> "Parameter "*string(targetdims[1])
    x = points_matrix[:,targetdims[1]]
    if sg.dims == 1
        seriestype  :=  :scatter
        x
    elseif  sg.dims == 2
        seriestype  :=  :scatter
        ylabel --> "Parameter "*string(targetdims[2])
        y = points_matrix[:,targetdims[2]]
        x,y
    elseif sg.dims > 2
        seriestype  :=  :scatter3d
        ylabel --> "Parameter "*string(targetdims[2])
        y = points_matrix[:,targetdims[2]]
        zlabel --> "Parameter "*string(targetdims[3])
        z = points_matrix[:,targetdims[3]]
        x,y,z
    end
end

"""
    f(sga::SparseGridApproximation)
    f(sg::SpectralSparseGridApproximation)

Recipe for plotting approximations based on sparse grids. For more than two dimensions it is a slice with the other parameters fixed at the midpoint.

"""
@recipe function f(sga::Union{SparseGridApproximation,SpectralSparseGridApproximation}; targetdims=[1,2])
    midpoint = 0.5*sum.(sga.domain)

    n = 100
    if length(midpoint) == 1
        subdomain1 = sga.domain[targetdims[1]]
        x =  range(subdomain1[1], stop=subdomain1[2], length=n)
        y = sga.(x)
        title --> "Sparse Grid Approximation"
        xlabel --> "Parameter "*string(targetdims[1])
        x,y
    else
        subdomain1 = sga.domain[targetdims[1]]
        subdomain2 = sga.domain[targetdims[2]]

        function  modifiedmidpoint(x,i,y,j)
            value = copy(midpoint)
            value[i] = x
            value[j] = y
            return value        
        end

        x =  range(subdomain1[1], stop=subdomain1[2], length=n)
        y =  range(subdomain2[1], stop=subdomain2[2], length=n)
        xy = [modifiedmidpoint(xi,targetdims[1],yi,targetdims[2]) for xi in x for yi in y]

        z = sga.(xy)

        seriestype  -->  :contour
        title --> "Sparse Grid Approximation"
        xlabel --> "Parameter "*string(targetdims[1])
        ylabel --> "Parameter "*string(targetdims[2])
        x,y,z
    end
end

