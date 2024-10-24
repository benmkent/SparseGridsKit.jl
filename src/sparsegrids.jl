using LinearAlgebra
using FastGaussQuadrature

import Base.:+
import Base.:-

mutable struct SparseGrid
    dims
    MI
    terms
    sparsegridptids
    maptermstosparsegrid
    cterms
    cmi
    pintsmap
    ptsperlevel
    ptsunique
    sparsegridpts
    polyperlevel
end

function Base.:+(sg1::SparseGrid, sg2::SparseGrid)
    dims = sg1.dims
    @assert dims == sg2.dims

    sg1maxmi = length(sg1.polyperlevel)
    sg2maxmi = length(sg2.polyperlevel)

    if sg1maxmi >= sg2maxmi
        sg = copy(sg1)
        sga = copy(sg1)
        sgb = copy(sg2)
    else
        sg = copy(sg2)
        sga = copy(sg2)
        sgb = copy(sg1)
    end

    for (imi, mi) in enumerate(eachcol(sgb.MI))
        index = findfirst(mi == mi2 for mi2 in eachcol(sga.MI))
        if index !== nothing
            sg.cmi[index] = sg.cmi[index] + sgb.cmi[imi]
        else
            push!(sg.cmi, sgb.cmi[imi])
            sg.MI = hcat(sg.MI, mi)
        end
    end

    # Clean zero c points
    #zeroindices = sg.cmi .!= 0
    #cmi = sg.cmi[zeroindices]
    #MI = sg.MI[:, zeroindices]

    (MI,cmi) = sortmi(sg.MI, sg.cmi)

    (terms, cterms, cmi, MI) = sparsegridterms(MI, sg.ptsperlevel, cmi=sg.cmi)

    (sparsegridptids, sparsegridpts, maptermstosparsegrid) = sparsegridpoints(terms, sg.pintsmap, sg.ptsunique, size(sg.MI, 1))

    sparsegrid = SparseGrid(dims, MI, terms, sparsegridptids, maptermstosparsegrid, cterms, cmi, sg.pintsmap, sg.ptsperlevel, sg.ptsunique, sparsegridpts, sg.polyperlevel)

    return sparsegrid
end

function Base.:-(sg1::SparseGrid, sg2::SparseGrid)
    sg2minus = copy(sg2)
    sg2minus.cterms = -sg2minus.cterms
    sg2minus.cmi = -sg2minus.cmi
    sg = copy(sg1) + sg2minus
    return sg
end

function Base.:copy(sg::SparseGrid)
    sgcopy = deepcopy(sg)
    return sgcopy
end

function sparsegridprecompute(maxmi)
    (ptsunique, ptsperlevel, polyperlevel, pintsmap) = sparsegridonedimdata(maxmi, doubling)
    productintegrals = sparsegridproductintegrals(polyperlevel, maxmi, doubling)
    return (; productintegrals)
end

function sparsegridproductintegrals(polyperlevel, maxmi, rule)
    @info "Computing product integrals"
    # Allocate arrays
    productintegrals = zeros(Float64, maxmi, maxmi, length(polyperlevel[end]), length(polyperlevel[end]))

    x, w = gausslegendre(2 * rule(maxmi))
    w = 0.5 * w
    wp1x = similar(x)

    pevalx = zeros(Float64, maxmi, length(polyperlevel[end]), length(x))
    for ii = 1:maxmi
        for jj = 1:length(polyperlevel[ii])
            pevalx[ii, jj, :] .= polyperlevel[ii][jj].(x)
        end
    end


    #productintegralsmatrix = Matrix{Float64}(undef, length(polyperlevel[end]), length(polyperlevel[end]))
    for ii = 1:maxmi
        # Compute integrals
        #productintegralslevel = []
        #pii = view(polyperlevel,ii);
        for ll = 1:ii
            #pll = view(polyperlevel,ll);s
            for p1 in eachindex(polyperlevel[ii])
                #p1x .= pii[1][p1].(x)
                #p1x .= w.* p1x
                wp1x = w .* pevalx[ii, p1, :]
                for p2 in eachindex(polyperlevel[ll])
                    #productintegralsmatrix[p1, p2] = sum(w .* polyperlevel[ii][p1].(x) .* polyperlevel[ll][p2].(x))
                    @inbounds productintegrals[ii, ll, p1, p2] = dot(wp1x, pevalx[ll, p2, :]) #p1x .* pll[1][p2].(x)
                    #display([ii,ll,p1,p2,sum(w .* polyperlevel[ii][p1].(x) .* polyperlevel[ll][p2].(x))])
                end
            end
            #push!(productintegralslevel, productintegralsmatrix[1:length(polyperlevel[ii]), 1:length(polyperlevel[ll])])
        end
        #push!(productintegrals, productintegralslevel)
    end
    return productintegrals
end

function createenhancedsparsegrid(sg)
    sgMI = downwards_closed_set(sg.MI)
    rm = reducedmargin(sgMI)
    enhancedMI = downwards_closed_set(hcat(sgMI, rm))
    enhancedsg = createsparsegrid(enhancedMI)
    return enhancedsg
end

function mapfromto(sgfrom, sgto)
    (fromptids, ~, ~) = sparsegridpoints(sgfrom.terms, sgto.pintsmap, sgto.ptsunique, sgto.dims)

    maptermstosparsegrid = Vector{Integer}(undef, size(fromptids, 1))
    for ii in 1:size(fromptids, 1)
        maptermstosparsegrid[ii] = findfirst(fromptids[ii, :] == uniquerow for uniquerow in eachrow(sgto.sparsegridptids))
    end
    return maptermstosparsegrid
end


function margin(MI)
    ndim = size(MI, 1)
    margin = zeros(Int32, ndim, 0)
    for mi in eachcol(MI)
        for i = 1:ndim
            ei = zeros(Int32, ndim)
            ei[i] = 1
            miplus = mi .+ ei
            if isnothing(findfirst(miplus == mi for mi in eachcol(MI)))
                margin = hcat(margin, miplus)
            end
        end
    end
    uniqueindices = unique(i -> margin[:, i], 1:size(margin, 2))
    return sortmi(margin[:, uniqueindices])
end

function sortmi(MI, cmi)

    sorted = sortmi([MI;transpose(cmi)])
    MI_sorted = sorted[1:end-1,:]
    cmi_sorted = copy(cmi)
    cmi_sorted .= sorted[end,:]
    return (MI_sorted, cmi_sorted)
end

function sortmi(MI)
    if ndims(MI) == 1
        sorted = sort(MI)
    else
        sorted = sortslices(MI, dims=2)
    end
    return sorted
end

function reducedmargin(MI)
    mgn = margin(MI)
    ndim = size(mgn, 1)
    rm = nothing

    for mi in eachcol(mgn)
        for ei in eachcol(I(ndim))
            miminus = mi .- ei
            if isnothing(findfirst(miminus == mi for mi in eachcol(MI)))
                break
            end
        end
        if isnothing(rm)
            rm = mi
        else
            rm = hcat(rm, mi)
        end
    end
    return sortmi(rm)
end

function createsparsegrid(MI; rule=doubling)
    MI = sortmi(MI)
    maxmi = maximum(MI)
    dims = size(MI, 1)

    (ptsunique, ptsperlevel, polyperlevel, pintsmap) = sparsegridonedimdata(maxmi, rule)

    (grid, cterms, cmi, MI) = sparsegridterms(MI, ptsperlevel)

    (sparsegridptids, sparsegridpts, maptermstosparsegrid) = sparsegridpoints(grid, pintsmap, ptsunique, dims)

    sparsegrid = SparseGrid(dims, MI, grid, sparsegridptids, maptermstosparsegrid, cterms, cmi, pintsmap, ptsperlevel, ptsunique, sparsegridpts, polyperlevel)
    return sparsegrid
end

function sparsegridpoints(grid, pintsmap, ptsunique, dims)
    gridasptsindices = Matrix{Integer}(undef, length(grid), dims)
    for ii = eachindex(grid)
        for jj = eachindex(grid[ii])
            gridasptsindices[ii, jj] = pintsmap[grid[ii][jj]]
        end
    end

    # Find unique points
    uniquerows = unique(i -> gridasptsindices[i, :], 1:size(gridasptsindices, 1))
    sparsegridptids = gridasptsindices[uniquerows, :]

    maptermstosparsegrid = Vector{Integer}(undef, size(gridasptsindices, 1))
    for ii in 1:size(gridasptsindices, 1)
        maptermstosparsegrid[ii] = findfirst(gridasptsindices[ii, :] == uniquerow for uniquerow in eachrow(sparsegridptids))
    end

    sparsegridpts = similar(sparsegridptids, Float64)
    for ii = 1:size(sparsegridptids, 1)
        for jj = 1:size(sparsegridptids, 2)
            sparsegridpts[ii, jj] = ptsunique[sparsegridptids[ii, jj]]
        end
    end
    return (sparsegridptids, sparsegridpts, maptermstosparsegrid)
end

function sparsegridterms(MI, ptsperlevel; cmi=nothing)
    if isnothing(cmi)
        cmi = computesparsegridc(MI)
    end
    terms = []
    cterms = []
    for (miind, mi) in enumerate(eachcol(MI))
        if true #cmi[miind] != 0
            indexarray = []
            for i in mi
                indexarrayi = []
                for j = eachindex(ptsperlevel[i])
                    push!(indexarrayi, [i, j])
                end
                push!(indexarray, collect(indexarrayi))
            end

            for row in collect(Iterators.product((indexarray[i] for i in eachindex(indexarray))...))
                push!(terms, collect(row))
                push!(cterms, cmi[miind])
            end
        end
    end
    # Clean up MI and cmi
    #zeroindices = cmi .!= 0
    #MI = MI[:, zeroindices]
    #cmi = cmi[zeroindices]

    return (terms, cterms, cmi, MI)
end

function sparsegridonedimdata(maxmi, rule)
    ptsunique = Vector{Float64}(undef, 0)
    ptsperlevel = Vector{Vector{Float64}}(undef, maxmi)
    #poly = Vector{Vector{Polynomial{Float64,:x}}}(undef, maxmi)
    polyperlevel = Vector{Any}(undef, maxmi)
    pintsmap = Dict{}()

    for ii = 1:maxmi
        ptsperlevel[ii] = ccpoints(rule(ii))
        polyperlevel[ii] = createlagrangepolys(ptsperlevel[ii])

        # Build points map
        for (jj, p) in enumerate(ptsperlevel[ii])
            if length(ptsunique) == 0
                ind = nothing
            else
                ind = findfirst(isapprox(p, ptsunique[i]) for i in eachindex(ptsunique))
            end

            if !isnothing(ind)
                pintsmap[[ii; jj]] = ind
            else
                push!(ptsunique, p)
                pintsmap[[ii; jj]] = length(ptsunique)
            end
        end
    end
    return (ptsunique, ptsperlevel, polyperlevel, pintsmap)
end

function interpolateonsparsegrid(sparsegrid, fongrid, targetpoints; vecdims=nothing)
    feval = Vector{Vector{ComplexF64}}(undef, size(targetpoints, 1))
    if isempty(fongrid)
        if isnothing(vecdims)
            error("Must input vector dimensions if using zero operator")
        end
        for k in eachindex(feval)
            feval[k] = zeros(vecdims)
        end
        return feval
    else
        vecdims = length(fongrid[1])
        for k in eachindex(feval)
            feval[k] = zeros(vecdims)
        end
    end
    terms = sparsegrid.terms
    cterms = sparsegrid.cterms
    maprowtouniquept = sparsegrid.maptermstosparsegrid
    poly = sparsegrid.polyperlevel

    val = [0.0]
    for k = 1:length(feval)
        targetpt = targetpoints[k, :]
        for (i, row) = enumerate(eachrow(terms))
            #val = c[i]*prod([polyperlevel[row[1][j][1]][row[1][j][2]](newpt[1][j]) for j in eachindex(row[1])])
            val .= cterms[i] * prod([
                poly[row[1][j][1]][row[1][j][2]](targetpt[j])
                for j in eachindex(row[1])])
            feval[k] = feval[k] .+ val[1] * fongrid[maprowtouniquept[i]]
        end
    end
    return feval
end

function integrateonsparsegrid(sparsegrid, fongrid, precompute; vecdims=nothing)
    if isempty(fongrid)
        if isnothing(vecdims)
            error("Must input vector dimensions if using zero operator")
        end
    else
        vecdims = length(fongrid[1])
    end
    terms = sparsegrid.terms
    cterms = sparsegrid.cterms
    maprowtouniquept = sparsegrid.maptermstosparsegrid
    productintegrals = precompute.productintegrals
    feval = Complex(0.0, 0.0) * ones(vecdims)
    val = [Complex(0.0, 0.0)]
    for (i, row) = enumerate(eachrow(terms))
        #val = c[i]*prod([polyperlevel[row[1][j][1]][row[1][j][2]](newpt[1][j]) for j in eachindex(row[1])])
        #val = cterms[i] * prod([productintegrals[row[1][j][1]][1][row[1][j][2], 1] for j in eachindex(row[1])])
        val .= cterms[i] .* prod(productintegrals[row[1][j][1], 1, row[1][j][2], 1] for j in eachindex(row[1]))
        feval .= feval .+ val .* fongrid[maprowtouniquept[i]]
    end
    return feval
end

function computepairwisenorms(fongrid, product)
    ngrid = length(fongrid)
    #pairwisenorms = Matrix(undef,ngrid, ngrid)
    pairwisenorms = fill(0.0 .* product(fongrid[1], fongrid[1]), (ngrid, ngrid))
    for i1 = 1:ngrid
        for i2 = 1:i1
            #pairwisenorms[i1, i2] = transpose(fongrid[i1]) * weight * fongrid[i2]
            pairwisenorms[i1, i2] = product(fongrid[i1], fongrid[i2])
        end
    end
    pairwisenorms .= pairwisenorms .+ conj(permutedims(pairwisenorms, [2, 1]))# .- Diagonal(diag(pairwisenorms))
    for ii = 1:ngrid
        pairwisenorms[ii, ii] = 0.5 .* pairwisenorms[ii, ii]
    end
    return pairwisenorms
end

function L2onsparsegrid(sparsegrid, fongrid, precompute; product=dot, pairwisenorms=nothing, returnpairwisenorms=false)
    if isnothing(pairwisenorms)
        @info "Must compute pairwise norms - may be possible to precompute for efficiency"
        pairwisenorms = computepairwisenorms(fongrid, product)
    end

    terms = sparsegrid.terms
    cterms = sparsegrid.cterms
    maprowtouniquept = sparsegrid.maptermstosparsegrid
    productintegrals = precompute.productintegrals
    feval = zeros(ComplexF64, size(pairwisenorms[1,1]))
    val = Complex(1.0, 0.0); #.* pairwisenorms[1, 1]]
    leveli = Vector{Integer}(undef, length(terms[1]))
    levelj = similar(leveli)
    indi = similar(leveli)
    indj = similar(leveli)

    valstomultiply = Vector{Float64}(undef,1)
    n = length(terms[1][1])
    for (i, rowi) = enumerate(eachrow(terms))
        if iszero(pairwisenorms[maprowtouniquept[i], :])
            continue
        end
        for (j, rowj) = enumerate(eachrow(terms))
            if iszero(pairwisenorms[maprowtouniquept[i], maprowtouniquept[j]])
                continue
            end
            @inbounds valstomultiply[1] = prod(productintegrals[rowi[1][k][1], rowj[1][k][1], rowi[1][k][2], rowj[1][k][2]] for k in eachindex(rowi[1]))
            # valstomultiply[1] = 1.0
            # for k in 1:n
            #     @inbounds a,b,c,d=rowi[1][k][1], rowj[1][k][1], rowi[1][k][2], rowj[1][k][2]
            #     @inbounds valstomultiply[1] *= productintegrals[a,b,c,d]
            # end
            @inbounds feval .= feval .+ cterms[i] * cterms[j] * valstomultiply[1] * pairwisenorms[maprowtouniquept[i], maprowtouniquept[j]]
        end
    end

    feval = sqrt.(abs.(feval)) # abs prevents very small negative numbers (≈ 0) giving error

    if returnpairwisenorms == true
        return (feval, pairwisenorms)
    else
        return feval
    end
end

function computesparsegridc(I)
    jjs = generatebinaryvectors(size(I, 1))
    cmi = Vector{Integer}(undef, size(I, 2))
    cmi .= 0
    iiplusjjs = similar(jjs)
    for (ii, i) = enumerate(eachcol(I))
        iiplusjjs = jjs .+ i
        for (jj, ij) = enumerate(eachcol(iiplusjjs))
            if any(ij == cmi for cmi in eachcol(I))
                cmi[ii] = cmi[ii] .+ (-1)^norm(jjs[:, jj], 1)
            end
        end
    end
    return cmi
end

function generatebinaryvectors(n)
    vectors = Matrix{Integer}(undef, n, 0)
    for i in 0:2^n-1
        binary_rep = reverse(digits(i, base=2, pad=n))
        vectors = hcat(vectors, binary_rep)
    end
    return vectors
end

function createlagrangepolys(pts)
    #p = Vector{Polynomial{Float64, :x}}(undef,length(pts))
    p = Vector{Any}(undef, length(pts))
    for ii in eachindex(pts)
        #p[ii] = fromroots(pts[1:end .!= ii])
        #p[ii] = p[ii]/p[ii](pts[ii])
        p[ii] = stablepoly(pts, ii)
    end
    return p
end

# function stablepoly(knots, ind)
#     n = length(knots)
#     l(z) = prod(z .- knots)
#     w_j(knots, j) = 1 / prod(knots[j] .- knots[1:n.!=j])

#     f_knots = zeros(length(knots))
#     f_knots[ind] = 1

#     w = [w_j(knots, ii) for ii in 1:n]

#     function L(x)
#         if any(x .== knots)
#             jj = findfirst(x .== knots)
#             y = f_knots[jj]
#         else
#             y = l(x) .* w[ind] ./ (x .- knots[ind])
#         end
#         return y
#     end
#     return L
# end

function stablepoly(knots, ind)
    n = length(knots)
    if n == 1
        function Lone(x)
            return 1.0
        end
        return Lone
    else
        l(z) = prod(z .- knots)

        f_knots = zeros(eltype(knots), length(knots))
        f_knots[ind] = 1

        function L(x)
            y = zero(eltype(knots))
            if any(x .== knots)
                jj = findfirst(x .== knots)
                y = f_knots[jj]
            else
                y = l(x) * 1 / prod(knots[ind] .- knots[ii] for ii in Iterators.filter(jj -> jj != ind, 1:n)) / (x .- knots[ind])
            end
            return y
        end
        return L
    end
end


doubling(n) = n == 1 ? 1 : 2^(n - 1) + 1

function createsmolyakmiset(n, k)
    iteratorI = fastvectorgenerator(k, n)
    I = Matrix{Int64}(undef, n, 0)
    for vector in iteratorI
        if norm(vector, 1) <= n + k
            I = hcat(I, collect(Int64,vector))
        end
    end
    I = sortmi(I)
    return I
end

fastvectorgenerator(k, n) = Iterators.product(repeat([1:(k+1)], n)...)

ccpoints(n) = n == 1 ? [0.0] : [0.5 * (cos(pi * (i / (n - 1)))) for i in 0:n-1] - [0.5 * (cos(pi * (i / (n - 1)))) for i in n-1:-1:0]

function is_leq(α::Vector, β::Vector)
    # Check if β ≤ α element-wise
    return all(β[i] ≤ α[i] for i in 1:length(α))
end

function downwards_closed_set(input_matrix::Matrix)
    num_indices, num_columns = size(input_matrix)
    output_matrix = copy(input_matrix)  # Initialize output matrix with input matrix

    # Iterate over each column (multi-index) in the input matrix
    for col_idx in 1:num_columns
        α = input_matrix[:, col_idx]

        # Generate all possible multi-indices β such that β ≤ α
        for β in Iterators.product((1:α[i] for i in 1:num_indices)...)
            β_vector = collect(β)
            if is_leq(α, β_vector) && !any(all(β_vector .== output_matrix[:, j]) for j in 1:size(output_matrix, 2))
                output_matrix = hcat(output_matrix, β_vector)
            end
        end
    end
    output_matrix = reduce(hcat, unique(eachcol(output_matrix)))

    return sortmi(output_matrix)
end
