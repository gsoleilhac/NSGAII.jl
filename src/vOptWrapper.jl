
function nsga(::Val{N}, popSize, nbGen, m, ϵ::Int = 5; kwargs...)::Vector{indiv{BitVector, Vector{Float64}, N, Float64}} where N
    vd = m.ext[:vOpt]
    @assert all(isfinite, m.colLower) "All variables must be bounded"
    @assert all(isfinite, m.colUpper) "All variables must be bounded"
    triplet = Tuple{Float64, Vector{Float64}, Vector{Int}}
    objs_triplets = SVector{N, triplet}([(obj.aff.constant, obj.aff.coeffs, map(x->getfield(x, :col), obj.aff.vars)) for obj in vd.objs])
    if all(equalto(:Min), vd.objSenses) || all(equalto(:Max), vd.objSenses)
        objSenses = SVector{N, Int}(fill(1, length(vd.objs)))
    else
        objSenses = SVector{N, Int}(map(s -> s == :Max ? -1 : 1, vd.objSenses)...)
    end
    _nsga(popSize, nbGen, m, vd, m.linconstr, objSenses, objs_triplets, ϵ ; kwargs...)
end

function nsga(popSize, nbGen, m, ϵ::Int = 5; kwargs...)
    vd = m.ext[:vOpt]
    @assert all(isfinite, m.colLower) "All variables must be bounded"
    @assert all(isfinite, m.colUpper) "All variables must be bounded"
    objs_triplets = SVector(((obj.aff.constant, obj.aff.coeffs, map(x->getfield(x, :col), obj.aff.vars)) for obj in vd.objs)...)
    if all(equalto(:Min), vd.objSenses) || all(equalto(:Max), vd.objSenses)
        objSenses = SVector{length(vd.objs)}(fill(1, length(vd.objs)))
    else
        objSenses = SVector(map(s -> s == :Max ? -1 : 1, vd.objSenses)...)
    end
    _nsga(popSize, nbGen, m, vd, m.linconstr, objSenses, objs_triplets, ϵ ; kwargs...)
end

function _nsga(popSize, nbGen, m, vd, linconstr, objSenses, objs_triplets::SVector{N, T}, ϵ = 5; kwargs...)::Vector{indiv{BitVector, Vector{Float64}, N, Float64}} where {N, T}

    mc = MixedCoding(ϵ, m.colCat, m.colLower, m.colUpper)
   
    function evaluate(cst, coeffs, vars, x)::Float64
        res = cst
        for i = 1:length(coeffs)
            @inbounds res += coeffs[i] * x[vars[i]]
        end
        res
    end

    z(x) = Tuple(objSenses .* map(obj_t->evaluate(obj_t..., x), objs_triplets))

    function CV(x)
        res::Float64 = 0.
        @inbounds for CSTR in linconstr

            if CSTR.lb != -Inf && CSTR.lb != typemin(Float64) 
                if CSTR.lb == 0
                    g = dot(CSTR.terms.coeffs, x[map(v-> getfield(v, :col), CSTR.terms.vars)])
                    res += max(0, -g)
                elseif CSTR.lb > 0 
                    g = dot(CSTR.terms.coeffs, x[map(v-> getfield(v, :col), CSTR.terms.vars)]) / CSTR.lb - 1
                    res += max(0, -g)
                else
                    g = dot(CSTR.terms.coeffs, x[map(v-> getfield(v, :col), CSTR.terms.vars)]) / CSTR.lb - 1
                    res += max(0, g)
                end

            end

            if CSTR.ub != Inf && CSTR.lb != typemax(Float64)
                 if CSTR.ub == 0
                    g = dot(CSTR.terms.coeffs, x[map(v-> getfield(v, :col), CSTR.terms.vars)])
                    res += max(0, g)
                elseif CSTR.ub > 0 
                    g = dot(CSTR.terms.coeffs, x[map(v-> getfield(v, :col), CSTR.terms.vars)]) / CSTR.ub - 1
                    res += max(0, g)
                else
                    g = dot(CSTR.terms.coeffs, x[map(v-> getfield(v, :col), CSTR.terms.vars)]) / CSTR.ub - 1
                    res += max(0, -g)
                end
            end
        end

        for i in find(m.colCat.==:Int)
            if x[i] > m.colUpper[i]
                res += m.colUpper - x[i]
            end
        end

        res
    end

    let mc=mc, vd=vd, z=z, CV=CV

        if all(equalto(:Max), vd.objSenses)
            res = nsga_max(popSize, nbGen, z, mc ; fCV = CV, kwargs...)
        else
            res = nsga(popSize, nbGen, z, mc ; fCV = CV, kwargs...)
        end

        if !all(equalto(:Min), vd.objSenses) || !all(equalto(:Max), vd.objSenses)
            for indiv in res
                indiv.y = indiv.y .* objSenses
            end
        end
        return res
    end
end
