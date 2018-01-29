module NSGAII
export nsga, MixedCoding, RealCoding
using ProgressMeter, Requires, StaticArrays

include("indivs.jl")
include("functions.jl")
include("crossover.jl")
include("mutation.jl")
include("mixedcoding.jl")


function nsga(popSize::Integer, nbGen::Integer, z::Function, init::Function ; 
    fdecode=identity, fCV = x->0., pmut= 0.05, fmut=default_mutation!, 
    fcross = default_crossover, seed=typeof(init())[], fplot = (x)->nothing)
    return _nsga(popSize, nbGen, init, z, fdecode, fCV , pmut, fmut, fcross, seed, fplot)
end

function nsga(popSize::Integer, nbGen::Integer, z::Function, mc::MixedCoding, init::Function=()->bitrand(mc.nbbitstotal); 
    fCV = x->0., pmut= 0.05, fmut=default_mutation!, fcross = default_crossover, 
    seed=Vector{Float64}[], fplot = (x)->nothing)
    return _nsga(popSize, nbGen, init, z, x->decode(x, mc), fCV , pmut, fmut, fcross, encode.(seed, mc), fplot)
end


function nsga(popSize, nbGen, m, ϵ = 5; kwargs...)
    vd = m.ext[:vOpt]
    @assert all(isfinite, m.colLower) "All variables must be bounded"
    @assert all(isfinite, m.colUpper) "All variables must be bounded"
    objs = SVector(vd.objs...)
    _nsga(popSize, nbGen, m, vd, objs, ϵ ; kwargs...)
end

function _nsga(popSize, nbGen, m, vd, objs::SVector{N, T}, ϵ = 5; kwargs...) where {N, T}

    mc = MixedCoding(ϵ, m.colCat, m.colLower, m.colUpper)
   
    evaluate(obj, x)::Float64 = evaluate(obj.aff.constant, obj.aff.coeffs, map(x->getfield(x, :col), obj.aff.vars), x)
    function evaluate(cst, coeffs, vars, x)::Float64
        res = cst
        for i = 1:length(coeffs)
            res += coeffs[i] * x[vars[i]]
        end
        res
    end
    z(x)::SVector{N, Float64} = map(obj->evaluate(obj, x), objs) .* SVector{N}(map(s -> s == :Max ? -1 : 1, vd.objSenses))

    function CV(x)
        res = 0.
        for CSTR in m.linconstr

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

    for i = 1:length(vd.objs)
        if vd.objSenses[i] == :Max
            vd.objs[i] = vd.objs[i] * -1
        end
    end

    let mc=mc, vd=vd, z=z, CV=CV

        if all(x->x==:Bin, m.colCat)
            res = nsga(popSize, nbGen, z, ()->bitrand(m.numCols) ; fCV = CV, kwargs...)
        else
            res = nsga(popSize, nbGen, z, mc ; fCV = CV, kwargs...)
        end

        for i = 1:length(vd.objs)
            if vd.objSenses[i] == :Max
                vd.objs[i] = vd.objs[i] * -1
            end
        end

        signs = SVector{N}([s == :Min ? 1 : -1 for s in vd.objSenses])

        for indiv in res
            indiv.y = indiv.y .* signs
        end

        res
        # [(ind.pheno, ind.y, ind.CV) for ind in res]

    end
end

end # module
