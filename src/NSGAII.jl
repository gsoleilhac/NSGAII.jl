module NSGAII

include("indivs.jl")
include("functions.jl")
include("crossover.jl")
include("mutation.jl")
include("realcoding.jl")

export nsga, RealCoding, decode, encode

using ProgressMeter

function nsga(popSize, nbGen, init, z, pMut= 0.05 ; fmut=default_mutation!, fcross = default_crossover, seed=typeof(init())[], fplot = (x)->nothing)
    P = [indiv(init(), z) for _=1:popSize]
    P[1:length(seed)] .= [indiv(s, z) for s in seed]
    fast_non_dominated_sort!(P)
    Q = similar(P)

    @showprogress 0.1 for gen = 1:nbGen
        
        for i = 1:popSize
            pa = tournament_selection(P)
            pb = tournament_selection(P)
            ca,cb = crossover(pa, pb, fcross)

            rand() < pMut && mutate!(ca, fmut)
            rand() < pMut && mutate!(cb, fmut)

            eval!(ca, z)
            eval!(cb, z)

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

end # module
