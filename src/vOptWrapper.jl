
function nsga(popSize, nbGen, m, ϵ::Int = 5; kwargs...)
    vd = m.ext[:vOpt]
    @assert all(isfinite, m.colLower) "All variables must be bounded"
    @assert all(isfinite, m.colUpper) "All variables must be bounded"
    objs_triplets = SVector(((obj.aff.constant, obj.aff.coeffs, map(x->getfield(x, :col), obj.aff.vars)) for obj in vd.objs)...)
    objSenses = SVector(map(s -> s == :Max ? -1 : 1, vd.objSenses)...)
    _nsga(popSize, nbGen, m, vd, m.linconstr, objSenses, objs_triplets, ϵ ; kwargs...)
end

function _nsga(popSize, nbGen, m, vd, linconstr, objSenses, objs_triplets::SVector{N, T}, ϵ = 5; kwargs...) where {N, T}

    mc = MixedCoding(ϵ, m.colCat, m.colLower, m.colUpper)
   
    function evaluate(cst, coeffs, vars, x)::Float64
        res = cst
        for i = 1:length(coeffs)
            @inbounds res += coeffs[i] * x[vars[i]]
        end
        res
    end
    z(x)::SVector{N, Float64} = objSenses .* map(obj_t->evaluate(obj_t..., x), objs_triplets)

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

        if all(x->x==:Bin, m.colCat)
            res = nsga(popSize, nbGen, z, ()->bitrand(m.numCols) ; fCV = CV, kwargs...)
        else
            res = nsga(popSize, nbGen, z, mc ; fCV = CV, kwargs...)
        end

        for indiv in res
            indiv.y = indiv.y .* objSenses
        end

        res

    end
end
