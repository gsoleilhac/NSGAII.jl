module NSGAII
export nsga, MixedCoding, RealCoding, decode, encode
using ProgressMeter, Requires

include("indivs.jl")
include("functions.jl")
include("crossover.jl")
include("mutation.jl")
include("mixedcoding.jl")


function nsga(popSize::Integer, nbGen::Integer, z::Function, init::Function ; fdecode=identity, fCV = x->0., pmut= 0.05, fmut=default_mutation!, fcross = default_crossover, seed=typeof(init())[], fplot = (x)->nothing)
    _nsga(popSize, nbGen, init, z, fdecode, fCV , pmut, fmut, fcross, seed, fplot)
end

function nsga(popSize::Integer, nbGen::Integer, z::Function, mc::MixedCoding, init::Function=()->bitrand(mc.nbbitstotal); fCV = x->0., pmut= 0.05, fmut=default_mutation!, fcross = default_crossover, seed=Vector{Float64}[], fplot = (x)->nothing)
    _nsga(popSize, nbGen, init, z, x->decode(x, mc), fCV , pmut, fmut, fcross, encode.(seed, mc), fplot)
end


@require vOptGeneric begin
function nsga(popSize, nbGen, m, ϵ = 5; kwargs...)

    vd = m.ext[:vOpt]
    @assert all(isfinite, m.colLower) "All variables must be bounded"
    @assert all(isfinite, m.colUpper) "All variables must be bounded"
    
    mc = MixedCoding(ϵ, m.colCat, m.colLower, m.colUpper)

    evaluate(obj, x) = dot(obj.aff.coeffs, x[map(v-> getfield(v, :col), obj.aff.vars)]) + obj.aff.constant

    z(x) = ((evaluate(obj, x) for obj in vd.objs)...)

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
        res = nsga(popSize, nbGen, z, mc ; fCV = CV, kwargs...)

        for i = 1:length(vd.objs)
            if vd.objSenses[i] == :Max
                vd.objs[i] = vd.objs[i] * -1
            end
        end

        signs = Tuple(s == :Min ? 1 : -1 for s in vd.objSenses)

        for indiv in res
            indiv.y = indiv.y .* signs
        end

        [(ind.pheno, ind.y, ind.CV) for ind in res]

    end
end
end

end # module
