using LinearAlgebra
using StaticArrays

import Base.:+
import Base.:-

mutable struct SparseGrid
    dims # Number of input dimensions
    MI  # Multi index set
    terms #
    sparsegridptids
    maptermstosparsegrid
    cterms
    cmi
    pintsmap
    ptsperlevel
    ptsunique
    sparsegridpts
    polyperlevel
    domain
end

mutable struct SparseGridApproximation
    sg::SparseGrid
    fongrid::Vector
end

function (sga::SparseGridApproximation)(x)
    @assert length(x) == sga.sg.dims
    return interpolate_on_sparsegrid(sga.sg, sga.fongrid, [x])[1]
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

    sparsegrid = SparseGrid(dims, MI, terms, sparsegridptids, maptermstosparsegrid, cterms, cmi, sg.pintsmap, sg.ptsperlevel, sg.ptsunique, sparsegridpts, sg.polyperlevel, sg.domain)

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

function sparsegridprecompute(maxmi, domain, knots=ccpoints, rule=doubling)
    if isa(knots, Function)
        knots = fill(knots, length(domain))
    end
    if isa(rule, Function)
        rule = fill(rule, length(domain))
    end
    productintegrals = Vector{Any}(undef, length(knots))
    for ii = eachindex(knots)
        (ptsunique, ptsperlevel, polyperlevel, pintsmap) = sparsegridonedimdata(maxmi, knots[ii], rule[ii], domain[ii])
        productintegrals[ii] = sparsegridproductintegrals(polyperlevel, maxmi, knots[ii], rule[ii])
    end
    @debug "Computed weighted L^2 product integrals in each dimension"
    return (; productintegrals)
end

function sparsegridproductintegrals(polyperlevel, maxmi, knots, rule)
    @debug "Computing product integrals"
    # Allocate arrays
    productintegrals = zeros(Float64, maxmi, maxmi, length(polyperlevel[end]), length(polyperlevel[end]))

    #x, w = gausslegendre(2 * rule(maxmi))
    x,w = knots(2*rule(maxmi))
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
        for ll = 1:maxmi
            #pll = view(polyperlevel,ll);s
            @inbounds for p1 in eachindex(polyperlevel[ii])
                #p1x .= pii[1][p1].(x)
                #p1x .= w.* p1x
                wp1x .= w .* pevalx[ii, p1, :]
                @inbounds for p2 in eachindex(polyperlevel[ll])
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

function mapfromto(sgfrom, sgto)
    (fromptids, ~, ~) = sparsegridpoints(sgfrom.terms, sgto.pintsmap, sgto.ptsunique, sgto.dims)

    maptermstosparsegrid = Vector(undef, size(fromptids, 1))
    for ii in 1:size(fromptids, 1)
        maptermstosparsegrid[ii] = findfirst(fromptids[ii, :] == uniquerow for uniquerow in eachrow(sgto.sparsegridptids))
    end
    @debug "Mapped each sparse grid point in sgfrom to sgto"
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
    rm = Matrix{Int64}(undef,ndim,0)

    for mi in eachcol(mgn)
        in_rm = true

        for ei in eachcol(I(ndim))
            miminus = mi .- ei
            
            test = any(v -> v == miminus, eachcol(MI)) || any(miminus .== 0)

            if test == false
                in_rm = false;
                break
            end
        end

        if in_rm == true
            rm = hcat(rm, mi)
        end
    end
    return sortmi(rm)
end

function createsparsegrid(MI, domain; rule=doubling, knots=ccpoints)
    MI = sortmi(MI)
    maxmi = maximum(MI)
    dims = size(MI, 1)

    @debug "Creating sparse grid with parameter dimension "*string(dims)

    if isa(knots, Function)
        knots = fill(knots, dims)
        @debug "Using the same notes in all dimensions"
    end
    if isa(rule, Function)
        rule = fill(rule, dims)
        @debug "Using the same rule in all dimensions"
    end
    try
        @assert (size(rule,1) == size(knots,1))
        @assert (size(rule,1) == dims)
    catch e
        DimensionMismatch("If specified, rule and knots must be the same length as the number of dimensions")
    end
    ptsunique = Vector{Any}(undef,dims)
    ptsperlevel = Vector{Any}(undef,dims)
    polyperlevel = Vector{Any}(undef,dims)
    pintsmap = Vector{Any}(undef,dims)
    for ii = eachindex(knots)
        (ptsunique[ii], ptsperlevel[ii], polyperlevel[ii], pintsmap[ii]) = sparsegridonedimdata(maxmi, knots[ii], rule[ii], domain[ii])
    end
    @debug "Constructed one dimensional grid data"

    (grid, cterms, cmi, MI) = sparsegridterms(MI, ptsperlevel)
    @debug "Constructed combintion technique data"

    filter_vector = cterms .!= 0
    grid = grid[filter_vector]
    cterms = cterms[filter_vector]

    (sparsegridptids, sparsegridpts, maptermstosparsegrid) = sparsegridpoints(grid, pintsmap, ptsunique, dims)
    @debug "Created sparse grid term to sparse grid points mappings"

    sparsegrid = SparseGrid(dims, MI, grid, sparsegridptids, maptermstosparsegrid, cterms, cmi, pintsmap, ptsperlevel, ptsunique, sparsegridpts, polyperlevel, domain)
    @debug "Created sparse grid"
    return sparsegrid
end

function sparsegridpoints(grid, pintsmap, ptsunique, dims)
    if ~isa(pintsmap,Vector)
        commonknots = true
    else
        commonknots = false
    end

    gridasptsindices = Matrix{Integer}(undef, length(grid), dims)
    for ii = eachindex(grid)
        for jj = eachindex(grid[ii])
            if commonknots == true
                gridasptsindices[ii, jj] = pintsmap[grid[ii][jj]]
            else
                gridasptsindices[ii, jj] = pintsmap[jj][grid[ii][jj]]
            end
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
            if commonknots == true
                sparsegridpts[ii, jj] = ptsunique[sparsegridptids[ii, jj]]
            else
                sparsegridpts[ii, jj] = ptsunique[jj][sparsegridptids[ii, jj]]
            end
        end
    end
    return (sparsegridptids, sparsegridpts, maptermstosparsegrid)
end

function sparsegridterms(MI, ptsperlevel; cmi=nothing)
    if isnothing(cmi)
        cmi = computesparsegridc(MI)
    end
    n = size(MI,1)
    terms = []
    cterms = []
    for (miind, mi) in enumerate(eachcol(MI))
        if true #cmi[miind] != 0
            indexarray = []
            for (ind,i) in enumerate(mi)
                indexarrayi = []
                if isa(ptsperlevel[1][1],Number)
                    for j = eachindex(ptsperlevel[i])
                        push!(indexarrayi, [i, j])
                    end
                else
                    for j = eachindex(ptsperlevel[ind][i])
                        push!(indexarrayi, [i, j])
                    end   
                end
                push!(indexarray, collect(indexarrayi))
            end

            for row in collect(Iterators.product((indexarray[i] for i in eachindex(indexarray))...))
                push!(terms,SVector{n}(row))
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

function sparsegridonedimdata(maxmi, knots, rule, domain)
    ptsunique = Vector{Real}(undef, 0)
    ptsperlevel = Vector{Vector{Real}}(undef, maxmi)
    #poly = Vector{Vector{Polynomial{Float64,:x}}}(undef, maxmi)
    polyperlevel = Vector{Any}(undef, maxmi)
    pintsmap = Dict{}()

    # Based on knots define appropriate ApproxFun space for underlying polynomials
    space = spaceselector(domain)

    for ii = 1:maxmi
        # Detect special case of fidelity
        if knots === fidelitypoints
            ptsperlevel[ii], w = knots(ii)
        else
            ptsperlevel[ii], w = knots(rule(ii))
        end
        polyperlevel[ii] = createlagrangepolys(ptsperlevel[ii], space)

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

function interpolateonsparsegrid(sparsegrid, fongrid, targetpoints; evaltype=nothing)
    if isempty(fongrid)
        if isnothing(evaltype)
            error("Must input evaltype if using zero operator")
        end
        feval = Vector{evaltype}(undef,length(targetpoints))
        for k in eachindex(feval)
            feval[k] .= zero(evaltype)
        end
        return feval
    else
        evaltype = typeof(fongrid[1])
        feval = Vector{evaltype}(undef,length(targetpoints))
    end
    terms = sparsegrid.terms
    cterms = sparsegrid.cterms
    maprowtouniquept = sparsegrid.maptermstosparsegrid
    # poly = sparsegrid.polyperlevel
    ndims = sparsegrid.dims
    if isa(sparsegrid.polyperlevel[1][1],Function)
        poly = Vector{Any}(undef,ndims) 
        for ii = 1:ndims
            poly[ii] = sparsegrid.polyperlevel
        end
    else
        poly = sparsegrid.polyperlevel 
    end


    cweightedpolyprod = zeros(length(feval),size(terms,1))
    for k = 1:length(feval), (i, row) = enumerate(eachrow(terms))
        targetpt = targetpoints[k]
        #for (i, row) = enumerate(eachrow(terms))
            #val = c[i]*prod([polyperlevel[row[1][j][1]][row[1][j][2]](newpt[1][j]) for j in eachindex(row[1])])
            try
            @inbounds  cweightedpolyprod[k,i] = cterms[i] * prod(poly[j][row[1][j][1]][row[1][j][2]](targetpt[j]) for j in eachindex(row[1]))
            #feval[k] = feval[k] .+ cweightedpolyprod[k,i] * fongrid[maprowtouniquept[i]]
            catch
                @show k, i
                @show cterms[i]
                @show row
                @show poly
                @show targetpt

                @show [poly[j][row[1][j][1]][row[1][j][2]](targetpt[j]) for j in eachindex(row[1])]
                rethrow()
            end
        #end
    end
    @inbounds @fastmath feval = [sum(cweightedpolyprod[k,i] * fongrid[maprowtouniquept[i]] for i in 1:size(terms,1)) for k=1:length(feval)]

    return feval
end

function integrateonsparsegrid(sparsegrid, fongrid, precompute; evaltype=nothing)
    if isempty(fongrid)
        if isnothing(evaltype)
            error("Must input evaltype if using zero operator")
        end
        fongrid[1] = 0.0*zero(evaltype)
    else
        evaltype = typeof(fongrid[1])
    end
    terms = sparsegrid.terms
    cterms = sparsegrid.cterms
    maprowtouniquept = sparsegrid.maptermstosparsegrid
    nterms = size(terms,1)
    ndims = sparsegrid.dims
    if ~isa(precompute.productintegrals,Vector)
        productintegrals = Vector{Any}(undef,ndims) 
        for ii = 1:ndims
            productintegrals[ii] = precompute.productintegrals
        end
    else
        productintegrals = precompute.productintegrals 
    end

    val = Vector{Float64}(undef,nterms)
    for (i, row) = enumerate(eachrow(terms))
        val[i] = cterms[i] * prod(productintegrals[j][row[1][j][1], 1, row[1][j][2], 1] for j in eachindex(row[1]))
    end
    @debug "Computed integral on sparse grid"
    feval = sum(val[i] * fongrid[maprowtouniquept[i]] for i=1:nterms)
    return feval
end

function computepairwisenorms(fongrid; product=nothing)
    if isnothing(product)
        product = (x,y) -> x.*y
    end
    ngrid = length(fongrid)
    #pairwisenorms = Matrix(undef,ngrid, ngrid)
    pairwisenorms = Matrix{typeof(product(fongrid[1], fongrid[1]))}(undef, (ngrid, ngrid))
    for i1 = 1:ngrid, i2 = 1:ngrid
        #for 
            #pairwisenorms[i1, i2] = transpose(fongrid[i1]) * weight * fongrid[i2]
            pairwisenorms[i1, i2] = product(fongrid[i1], fongrid[i2])
        #end
    end
    # pairwisenorms .= pairwisenorms .+ conj(permutedims(pairwisenorms, [2, 1]))# .- Diagonal(diag(pairwisenorms))
    # for ii = 1:ngrid
    #     pairwisenorms[ii, ii] = 0.5 .* pairwisenorms[ii, ii]
    # end
    return pairwisenorms
end

function L2onsparsegrid(sparsegrid, fongrid, precompute; product=nothing, pairwisenorms=nothing, returnpairwisenorms=false)
    if isnothing(pairwisenorms)
        @debug "Must compute pairwise norms - may be possible to precompute for efficiency"
        pairwisenorms = computepairwisenorms(fongrid, product=product)
    end

    terms = sparsegrid.terms
    cterms = sparsegrid.cterms
    maprowtouniquept = sparsegrid.maptermstosparsegrid
    ndims = sparsegrid.dims
    if ~isa(precompute.productintegrals,Vector)
        productintegrals = Vector{Any}(undef,ndims) 
        for ii = 1:ndims
            productintegrals[ii] = precompute.productintegrals
        end
    else
        productintegrals = precompute.productintegrals 
    end
    feval = zeros(ComplexF64, size(pairwisenorms[1,1]))

    nonzero_rows = [any(!iszero(pairwisenorms[maprowtouniquept[i], :])) for i in 1:size(terms,1)]
    valstomultiply = ones(size(terms,1),size(terms,1))
    n = length(terms[1][1])

    a = Vector{Int64}(undef,2)
    b = Vector{Int64}(undef,2)

    for (i, rowi) = enumerate(eachrow(terms)), (j, rowj) = enumerate(eachrow(terms))
        @inbounds if nonzero_rows[i] == false
            continue
        end
        #for (j, rowj) = enumerate(eachrow(terms))
            #@inbounds if iszero(pairwisenorms[maprowtouniquept[i], maprowtouniquept[j]])
            #    display(nonzero_rows[j])
            @inbounds if nonzero_rows[j] == false
                continue
            end
            #@inbounds @fastmath valstomultiply[i,j] = prod(productintegrals[rowi[1][k][1], rowj[1][k][1], rowi[1][k][2], rowj[1][k][2]] for k in eachindex(rowi[1]))
            # valstomultiply[1] = 1.0
            rowi1 = rowi[1]
            rowj1 = rowj[1]
            for k in eachindex(rowi1)
                @inbounds a .= rowi1[k]
                @inbounds b .= rowj1[k]
                @inbounds @fastmath valstomultiply[i,j] *= productintegrals[k][a[1],b[1],a[2],b[2]]
            end
            feval .= feval .+ cterms[i] * cterms[j] * valstomultiply[i,j] * pairwisenorms[maprowtouniquept[i], maprowtouniquept[j]]
        #end
    end
    #@inbounds feval .= sum(sum(cterms[i] * cterms[j] * valstomultiply[i,j] * pairwisenorms[maprowtouniquept[i], maprowtouniquept[j]] for j=1:size(terms,1)) for i = 1:size(terms,1))

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

function createlagrangepolys(pts, space=Chebyshev(-1.0..1.0))
    p = Vector{Any}(undef, length(pts))
    for ii in eachindex(pts)
        # p[ii] = stablepoly(pts, ii)
        p[ii] = poly_Fun(pts, ii, space=space)
    end
    return p
end

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

function createtensormiset(n, k)
    iteratorI = fastvectorgenerator(k, n)
    I = Matrix{Int64}(undef, n, 0)
    for vector in iteratorI
        if maximum(vector) <= n + k
            I = hcat(I, collect(Int64,vector))
        end
    end
    I = sortmi(I)
    return I
end

fastvectorgenerator(k, n) = Iterators.product(repeat([1:(k+1)], n)...)

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

function spaceselector(domain)
    if length(domain) != 2
        Base.throw_checksize_error("Domain must be 2 values")
    end
    xmin = domain[1]
    xmax = domain[2]
    domain = xmin..xmax

    if isinf(domain.left) && isinf(domain.right)
        space = Hermite(domain)
    elseif iszero(domain.left) && isinf(domain.right)
        space = Laguerre()
    else
        space = Chebyshev(domain)
    end
    return space
end
