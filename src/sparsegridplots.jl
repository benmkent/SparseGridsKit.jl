using RecipesBase
using StatsBase

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
    misetplot(miset::MISet)

    Scatter plot of MI in pairs of dimensions

# Arguments
- `miset` : Multi-index set to plot
"""
@userplot MISetPlot
@recipe function f(h::MISetPlot)
    miset = h.args[1]
    @assert isa(miset, MISet)
    data = hcat(get_mi(miset)...)'

    n, k = size(data, 1), size(data, 2)
    
    title --> "Multi-Index Pair Plot"

    # Set up subplots (creating a k x k grid)
    layout := @layout (k,k)  # Create a k x k layout grid

    # Loop through the data to create scatter plots
    for jj = 1:k
        for ii = 1:jj  # Only plot the lower triangle and diagonal
            x, y = data[:, ii], data[:, jj]
            # # Marker size for frequency
            # points = [(x[i], y[i]) for i in 1:length(x)]
            # point_counts = countmap(points)
            # unique_points = keys(point_counts)
            # frequencies = values(point_counts)
            # sizes = [100 * frequencies[findfirst(isequal(p), unique_points)] for p in points]
            @series begin
                subplot := ii + k*(jj-1)  # Assign subplot to the correct position
                seriestype := :scatter  # Define the plot type
                xlabel --> "Variable $ii"
                ylabel --> "Variable $jj"
                # Return x,y
                x,y
            end
        end
    end
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
    x,y,z = nothing, nothing, nothing

    seriestype  :=  :scatter3d
    title --> "Sparse Grid"
    xlabel --> "Parameter "*string(targetdims[1])
    x = points_matrix[:,targetdims[1]]
    y = zeros(size(x))
    z = zeros(size(x))
    if sg.dims > 1
        ylabel --> "Parameter "*string(targetdims[2])
        y = points_matrix[:,targetdims[2]]
    end
    if sg.dims > 2
        zlabel --> "Parameter "*string(targetdims[3])
        z = points_matrix[:,targetdims[3]]
    end
    x,y,z
end

"""
    f(sga::SparseGridApproximation)
    f(sg::SpectralSparseGridApproximation)

Recipe for plotting approximations based on sparse grids. For more than two dimensions it is a slice with the other parameters fixed at the midpoint.

"""
@recipe function f(sga::Union{SparseGridApproximation,SpectralSparseGridApproximation}; targetdims=[1,2])
    if sga.sg.dims == 1
        title --> "Sparse Grid Approximation"
        xlabel --> "Parameter "*string(targetdims[1])

        subdomain1 = sga.domain[targetdims[1]]

        n = 100
        infreplacement = 10
        a1 = isinf(subdomain1[1]) ? -infreplacement : subdomain1[1]
        b1 = isinf(subdomain1[2]) ? infreplacement : subdomain1[2]
        x = range(a1, stop=b1, length=n)
        y = sga.(x)
        
        x,y
    else
        midpoint = mean.(sga.domain)

        subdomain1 = sga.domain[targetdims[1]]
        subdomain2 = sga.domain[targetdims[2]]

        function  modifiedmidpoint(x,i,y,j)
            value = copy(midpoint)
            value[i] = x
            value[j] = y
            return value        
        end

        n = 2
        infreplacement = 10
        a1 = isinf(subdomain1[1]) ? -infreplacement : subdomain1[1]
        b1 = isinf(subdomain1[2]) ? infreplacement : subdomain1[2]
        a2 = isinf(subdomain2[1]) ? -infreplacement : subdomain2[1]
        b2 = isinf(subdomain2[2]) ? infreplacement : subdomain2[2]
        x = range(a1, stop=b1, length=n)
        y = range(a2, stop=b2, length=n)
        xy = [modifiedmidpoint(xi,targetdims[1],yi,targetdims[2]) for yi in y for xi in x]

        print(collect(xy))

        z = sga.(xy)

        print(z)

        seriestype  -->  :contour
        title --> "Sparse Grid Approximation"
        xlabel --> "Parameter "*string(targetdims[1])
        ylabel --> "Parameter "*string(targetdims[2])

        x,y,z
    end
end

