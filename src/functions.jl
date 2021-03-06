function _nsga(::Indiv{G, Y}, sense, popSize, nbGen, init, z, 
    fCV, pmut, fmut, fcross, seed, fplot, plotevery, refreshtime)::Vector{Indiv{G, Y}} where {G,Y}

    popSize = max(popSize, length(seed))
    isodd(popSize) && (popSize += 1)
    P = Vector{Indiv{G, Y}}(undef, 2 * popSize)
    P[1:(popSize - length(seed))] .= [createIndiv(init(), z, fCV) for _ = 1:(popSize - length(seed))]
    for i = 1:length(seed)
        P[popSize - length(seed) + i] = createIndiv(convert(G, seed[i]), z, fCV)
        if fCV(P[popSize - length(seed) + i].x) > 0
            @warn "element $i of the seed is unfeasible"
        end
    end
    for i = 1:popSize
        P[popSize+i] = deepcopy(P[i])
    end
    fast_non_dominated_sort!(view(P, 1:popSize), sense)

    @showprogress refreshtime for gen = 1:nbGen
        for i = 1:2:popSize

            pa = tournament_selection(P)
            pb = tournament_selection(P)

            crossover!(pa, pb, fcross, P[popSize + i], P[popSize + i + 1])

            rand() < pmut && mutate!(P[popSize + i], fmut)
            rand() < pmut && mutate!(P[popSize + i + 1], fmut)

            eval!(P[popSize + i], z, fCV)
            eval!(P[popSize + i + 1], z, fCV)
        end

        fast_non_dominated_sort!(P, sense)
        sort!(P, by = x -> x.rank, alg = Base.Sort.QuickSort)
        
        let f::Int = 1
            ind = 0
            indnext = findlast(x -> x.rank == f, P)
            while 0 < indnext <= popSize
                ind = indnext
                f += 1
                indnext = findlast(x -> x.rank == f, P)
            end
            indnext == 0 && (indnext = length(P))
            crowding_distance_assignment!(view(P, ind+1:indnext))
            sort!(view(P, (ind + 1):indnext), by = x -> x.crowding, rev = true, alg = PartialQuickSort(popSize - ind))
        end

        gen % plotevery == 0 && fplot(P)
    end
    fplot(P)
    filter(x -> x.rank == 1, view(P, 1:popSize))
end

function fast_non_dominated_sort!(pop::AbstractVector{T}, sense) where {T}
    n = length(pop)

    for p in pop
        empty!(p.dom_list)
        p.dom_count = 0
        p.rank = 0
    end

    @inbounds for i in 1:n
        for j in i+1:n
            if dominates(sense, pop[i], pop[j])
                push!(pop[i].dom_list, j)
                pop[j].dom_count += 1
            elseif dominates(sense, pop[j], pop[i])
                push!(pop[j].dom_list, i)
                pop[i].dom_count += 1
            end
        end
        if pop[i].dom_count == 0
            pop[i].rank = 1
        end
    end

    k = UInt16(2)
    @inbounds while any(==(k-one(UInt16)), (p.rank for p in pop)) #ugly workaround for #15276
        for p in pop 
            if p.rank == k-one(UInt16)
                for q in p.dom_list
                    pop[q].dom_count -= one(UInt16)
                    if pop[q].dom_count == zero(UInt16)
                        pop[q].rank = k
                    end
                end
            end
        end
        k += one(UInt16)
    end
    nothing
end

function crowding_distance_assignment!(pop::AbstractVector{Indiv{X, NTuple{N, T}}}) where {X, N, T}
    if N == 2
        sort!(pop, by = x -> x.y[1])
        pop[1].y[1] == pop[end].y[1] && return #Don't waste time if all indivs are the same
        pop[1].crowding = pop[end].crowding = Inf

        width_y1 = (pop[end].y[1] - pop[1].y[1])
        width_y2 = (pop[1].y[2] - pop[end].y[2])
        @inbounds for i = 2:length(pop)-1
            pop[i].crowding = (pop[i+1].y[1] - pop[i-1].y[1]) / width_y1 + (pop[i-1].y[2] - pop[i+1].y[2]) / width_y2
        end
    else
        for ind in pop
            ind.crowding = 0.
        end
        @inbounds for j = 1:length(first(pop).y) # Foreach objective
            let j = j #https://github.com/JuliaLang/julia/issues/15276
                sort!(pop, by = x -> x.y[j]) #sort by the objective value
            end
            pop[1].crowding = pop[end].crowding = Inf #Assign infinite value to extremas
            if pop[1].y[j] != pop[end].y[j]
                for i = 2:length(pop)-1
                    pop[i].crowding += (pop[i+1].y[j] - pop[i-1].y[j]) / (pop[end].y[j] - pop[1].y[j])
                end
            end
        end
    end
end

function tournament_selection(P)
    a, b = rand(1:length(P)÷2), rand(1:length(P)÷2)
    P[a] < P[b] ? P[a] : P[b]
end