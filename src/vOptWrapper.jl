function nsga(popSize, nbGen, m, ϵ::Int = 5; kwargs...)::Vector{indiv{BitVector, Vector{Float64}, Vector{Float64}}}
    vd = m.ext[:vOpt]
    @assert all(isfinite, m.colLower) "All variables must be bounded"
    @assert all(isfinite, m.colUpper) "All variables must be bounded"
    objs_constant = [obj.aff.constant for obj in vd.objs]
    objs_coeffs = [obj.aff.coeffs for obj in vd.objs]
    objs_vars = [map(x->getfield(x, :col), obj.aff.vars) for obj in vd.objs]
    if all(equalto(:Min), vd.objSenses) || all(equalto(:Max), vd.objSenses)
        objSenses = fill(1, length(vd.objs))
    else
        objSenses = map(s -> s == :Max ? -1 : 1, vd.objSenses)
    end
    _nsga(popSize, nbGen, m, vd, m.linconstr, objSenses, objs_constant, objs_coeffs, objs_vars, ϵ ; kwargs...)
end

function nsga_binary(popSize, nbGen, m ; kwargs...)::Vector{indiv{BitVector, BitVector, Vector{Float64}}}
    vd = m.ext[:vOpt]
    objs_constant = [obj.aff.constant for obj in vd.objs]
    objs_coeffs = [obj.aff.coeffs for obj in vd.objs]
    objs_vars = [map(x->getfield(x, :col), obj.aff.vars) for obj in vd.objs]
    if all(equalto(:Min), vd.objSenses) || all(equalto(:Max), vd.objSenses)
        objSenses = fill(1, length(vd.objs))
    else
        objSenses = map(s -> s == :Max ? -1 : 1, vd.objSenses)
    end
    _nsga_binary(popSize, nbGen, m, vd, m.linconstr, objSenses, objs_constant, objs_coeffs, objs_vars ; kwargs...)
end


function _nsga(popSize, nbGen, m, vd, linconstr, objSenses, objs_constant, objs_coeffs, objs_vars, ϵ = 5; kwargs...)::Vector{indiv{BitVector, Vector{Float64}, Vector{Float64}}}

    mc = MixedCoding(ϵ, m.colCat, m.colLower, m.colUpper)
   
    function evaluate(cst, coeffs, vars, x)::Float64
        res = cst
        for i = 1:length(coeffs)
            @inbounds res += coeffs[i] * x[vars[i]]
        end
        res
    end

    z(x) = objSenses .* map(i->evaluate(objs_constant[i], objs_coeffs[i], objs_vars[i], x), 1:length(objSenses))

    cstr_var_indices = [map(v-> getfield(v, :col), CSTR.terms.vars) for CSTR in linconstr]
    cstr_coeffs = [CSTR.terms.coeffs for CSTR in linconstr]
    cstr_lb = [CSTR.lb for CSTR in linconstr]
    cstr_ub = [CSTR.ub for CSTR in linconstr]
    intvars = find(m.colCat.==:Int)

    function CV(x)
        res::Float64 = 0.
        @inbounds for i in eachindex(linconstr)

            if cstr_lb[i] != -Inf && cstr_lb[i] != typemin(Float64) 
                if cstr_lb[i] == 0
                    g = dot(cstr_coeffs[i], view(x, cstr_var_indices[i]))
                    res += max(0, -g)
                elseif cstr_lb[i] > 0 
                    g = dot(cstr_coeffs[i], view(x, cstr_var_indices[i])) / cstr_lb[i] - 1
                    res += max(0, -g)
                else
                    g = dot(cstr_coeffs[i], view(x, cstr_var_indices[i])) / cstr_lb[i] - 1
                    res += max(0, g)
                end

            end

            if cstr_ub[i] != Inf && cstr_ub[i] != typemax(Float64)
                 if cstr_ub[i] == 0
                    g = dot(cstr_coeffs[i], view(x, cstr_var_indices[i]))
                    res += max(0, g)
                elseif cstr_ub[i] > 0 
                    g = dot(cstr_coeffs[i], view(x, cstr_var_indices[i])) / cstr_ub[i] - 1
                    res += max(0, g)
                else
                    g = dot(cstr_coeffs[i], view(x, cstr_var_indices[i])) / cstr_ub[i] - 1
                    res += max(0, -g)
                end
            end
        end

        for i in intvars
            if x[i] > m.colUpper[i]
                res += m.colUpper - x[i]
            end
        end

        res
    end

    let mc=mc

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

function _nsga_binary(popSize, nbGen, m, vd, linconstr, objSenses, objs_constant, objs_coeffs, objs_vars ; kwargs...)::Vector{indiv{BitVector, BitVector, Vector{Float64}}}
   
    function evaluate(cst, coeffs, vars, x)::Float64
        res = cst
        for i = 1:length(coeffs)
            @inbounds res += coeffs[i] * x[vars[i]]
        end
        res
    end

    z(x) = objSenses .* map(i->evaluate(objs_constant[i], objs_coeffs[i], objs_vars[i], x), 1:length(objSenses))

    cstr_var_indices = [map(v-> getfield(v, :col), CSTR.terms.vars) for CSTR in linconstr]
    cstr_coeffs = [CSTR.terms.coeffs for CSTR in linconstr]
    cstr_lb = [CSTR.lb for CSTR in linconstr]
    cstr_ub = [CSTR.ub for CSTR in linconstr]

    function CV(x)
        res::Float64 = 0.
        @inbounds for i in eachindex(linconstr)

            if cstr_lb[i] != -Inf && cstr_lb[i] != typemin(Float64) 
                if cstr_lb[i] == 0
                    g = dot(cstr_coeffs[i], view(x, cstr_var_indices[i]))
                    res += max(0, -g)
                elseif cstr_lb[i] > 0 
                    g = dot(cstr_coeffs[i], view(x, cstr_var_indices[i])) / cstr_lb[i] - 1
                    res += max(0, -g)
                else
                    g = dot(cstr_coeffs[i], view(x, cstr_var_indices[i])) / cstr_lb[i] - 1
                    res += max(0, g)
                end

            end

            if cstr_ub[i] != Inf && cstr_ub[i] != typemax(Float64)
                 if cstr_ub[i] == 0
                    g = dot(cstr_coeffs[i], view(x, cstr_var_indices[i]))
                    res += max(0, g)
                elseif cstr_ub[i] > 0 
                    g = dot(cstr_coeffs[i], view(x, cstr_var_indices[i])) / cstr_ub[i] - 1
                    res += max(0, g)
                else
                    g = dot(cstr_coeffs[i], view(x, cstr_var_indices[i])) / cstr_ub[i] - 1
                    res += max(0, -g)
                end
            end
        end

        res
    end

    # x = bitrand(m.numCols)
    # @code_warntype CV(x)

    let vd=vd

        if all(equalto(:Max), vd.objSenses)
            res = nsga_max(popSize, nbGen, z, ()->bitrand(m.numCols) ; fCV = CV, kwargs...)
        else
            res = nsga(popSize, nbGen, z, ()->bitrand(m.numCols) ; fCV = CV, kwargs...)
        end

        if !all(equalto(:Min), vd.objSenses) || !all(equalto(:Max), vd.objSenses)
            for indiv in res
                indiv.y = indiv.y .* objSenses
            end
        end
        return res
    end
end
