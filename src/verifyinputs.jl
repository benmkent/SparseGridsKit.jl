"""
    verifyinputs!(miset::MISet; forceadmissibility=true)

Checks miset for admissibility

# Arguments
- `miset`: Ensures miset is admissible.
- `forceadmissibility`: (optional) updates miset to be admissible

# Returns
- `miset`: Updated miset
"""
function verifyinputs!(miset::MISet; forceadmissibility=true)
    MI = get_mi(miset)

    n = length(MI[1])

    checkall = true
    
    # Needs fixing
    # The MI set is not necessarily complete. MI with combination techniuqe coefficients equal to 0 are removed.
    # for mi in MI
    #     backwardneighbours = mi .- I(n)
    #     for backmi in eachcol(backwardneighbours)
    #         backmi = max.(backmi,ones(Integer,n))
    #         check = any(backmi .== MI)

    #         if check == true
    #         elseif forceadmissibility
    #             miset = add_mi(miset, backmi)
    #             MI = get_mi(miset)
    #             @warn "Multi-index set is not admissible. Adding multi-index "*string(backmi)*" to ensure admissibility."
    #         else
    #             @error "Multi-index set is not admissible. Missing "*string(backmi)*"."
    #             checkall = false
    #         end
    #     end
    # end
    return checkall
end

"""
    verifyinputs(domain, knots)

Checks domain and knots are compatible

# Arguments
- `domain`: Approximation domain.
- `knots`: Knots function (or vector of knots functions)

# Returns
- `check`: Boolean
"""
function verifyinputs(domain, knots::Union{Function,Vector{Function}})
    testlevel = 3
    check = true
    testknots = Vector{Vector{Float64}}(undef, length(domain))
    for ii in eachindex(domain)
        if typeof(knots) == Vector{Function}
            (x,w) =  knots[ii](testlevel)
        else
            (x,w) = knots(testlevel)
        end
        testknots[ii] = [minimum(x), maximum(x)]
        if !isapprox(sum(w),1.0)
            @warn "Weights do not sum to 1."
        end
    end
    combinations = Iterators.product(testknots...)
    combinations = collect(combinations)
    for ii in eachindex(domain)
        if verifyinputs(domain, combinations) == false
            @error "Knots are not compatible with domain."
            check = false
        end
    end
    return check
end

"""
    verifyinputs(interval, points)

Checks points are within domain

# Arguments
- `domain`: Domain.
- `points`: Points to check.

# Returns
- `check`: Boolean
"""
function verifyinputs(domain, points)
    check = true
    for p in points
        if all(p[ii] < domain[ii][1] || p[ii] > domain[ii][2] for ii in eachindex(domain))
            @warn "Point "*string(p)*" is not within domain."
            check = false
            break
        end
    end
    return check
end