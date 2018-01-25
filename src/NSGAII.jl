module NSGAII

include("indivs.jl")
include("functions.jl")
include("crossover.jl")
include("mutation.jl")
include("realcoding.jl")

export nsga, RealCoding, decode, encode

using ProgressMeter, Requires

function nsga(popSize::Integer, nbGen::Integer, init, z ; fCV = x->0., pmut= 0.05, fmut=default_mutation!, fcross = default_crossover, seed=typeof(init())[], fplot = (x)->nothing)


    X = typeof(init())
    P = [indiv(init(), z, fCV) for _=1:popSize-length(seed)]
    append!(P, indiv.(convert.(X, seed),z, fCV))
    fast_non_dominated_sort!(P)
    Q = similar(P)

    @showprogress 0.1 for gen = 1:nbGen
        
        for i = 1:popSize
            pa = tournament_selection(P)
            pb = tournament_selection(P)
            ca,cb = crossover(pa, pb, fcross)

            rand() < pmut && mutate!(ca, fmut)
            rand() < pmut && mutate!(cb, fmut)

            eval!(ca, z, fCV)
            eval!(cb, z, fCV)

            if ca ⋖ cb
                Q[i] = ca
            elseif cb ⋖ ca
                Q[i] = cb
            else
                Q[i] = ifelse(rand(Bool), ca, cb)
            end
        end

        F = fast_non_dominated_sort!(vcat(P, Q))
        i = 1
        empty!(P)
        while length(P) + length(F[i]) <= popSize
            append!(P, F[i])
            i += 1
        end
        if length(P) != popSize
            crowding_distance_assignement!(F[i])
            sort!(F[i], by = x -> x.crowding)
            while length(P) < popSize
                push!(P, pop!(F[i]))
            end
        end

        fplot(P)
    end
    P
end


@require vOptGeneric begin
function nsga(popSize, nbGen, m ; kwargs...)

    vd = @eval Main getvOptData(m)
    @assert all(isfinite, m.colLower) "All variables must be bounded"
    @assert all(isfinite, m.colUpper) "All variables must be bounded"

    @assert !(:Int in m.colCat) "Only continuous and binary variables are supported"
    ϵ = map(x -> x==:Cont ? 5 : 0, m.colCat)
    d = RealCoding(ϵ, m.colLower, m.colUpper)

    init = () -> bitrand(d.nbbitstotal)

    function evaluate(obj, x)
        dot(obj.aff.coeffs, x[map(v-> getfield(v, :col), obj.aff.vars)]) + obj.aff.constant
    end

    function z(bits)
        x = decode(bits, d)
        ((evaluate(obj, x) for obj in vd.objs)...)
    end

    function CV(bits)
        x = decode(bits, d)
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
        res
    end

    for i = 1:length(vd.objs)
        if vd.objSenses[i] == :Max
            vd.objs[i] = vd.objs[i] * -1
        end
    end

    @code_warntype nsga(popSize, nbGen, init, z ; fCV = CV, kwargs...)
    res = nsga(popSize, nbGen, init, z ; fCV = CV, kwargs...)

    for i = 1:length(vd.objs)
        if vd.objSenses[i] == :Max
            vd.objs[i] = vd.objs[i] * -1
        end
    end

    signs = Tuple(s == :Min ? 1 : -1 for s in vd.objSenses)

    @show signs

    for indiv in res
        indiv.y = indiv.y .* signs
    end

    [(decode(ind.x, d), ind.y, ind.CV) for ind in res]

end
end

end # module
