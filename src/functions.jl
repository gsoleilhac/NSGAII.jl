function _nsga(::Type{indiv{G,Ph,N,Y}}, popSize, nbGen, init, z, fdecode, fCV , pmut, fmut, fcross, seed, fplot) where {G,Ph,N,Y}

    P::Vector{indiv{G,Ph,N,Y}} = Vector{indiv{G,Ph,N,Y}}(uninitialized, 2*popSize)
    P[1:popSize-length(seed)] .= [indiv(init(), fdecode, z, fCV) for _=1:popSize-length(seed)]
    for i = 1:length(seed)
        P[popSize-length(seed)+i] = indiv(convert(G, seed[i]), fdecode, z, fCV)
    end
    for i=1:popSize
        P[popSize+i] = deepcopy(P[i])
    end
    fast_non_dominated_sort!(view(P, 1:popSize))

    @showprogress 0.2 for gen = 1:nbGen
        
        for i = 1:2:popSize
            pa = tournament_selection(view(P, 1:popSize))
            pb = tournament_selection(view(P, 1:popSize))

            crossover!(pa, pb, fcross, P[popSize+i], P[popSize+i+1])

            rand() < pmut && mutate!(P[popSize+i], fmut)
            rand() < pmut && mutate!(P[popSize+i+1], fmut)

            eval!(P[popSize+i], fdecode, z, fCV)
            eval!(P[popSize+i+1], fdecode, z, fCV)
        end

        fast_non_dominated_sort!(P)
        sort!(P, by = x->x.rank, alg=Base.Sort.QuickSort)
        
        let f::Int = 1
            ind = 0
            indnext = findlast(x->x.rank==f, P)
            while 0 < indnext <= popSize
                ind = indnext
                f += 1
                indnext = findlast(x->x.rank==f, P)
            end
            indnext == 0 && (indnext = length(P))
            crowding_distance_assignement!(view(P, ind+1:indnext))
            sort!(view(P, ind+1:indnext), by = x -> x.crowding, rev=true, alg=PartialQuickSort(popSize-ind))
        end

        fplot(P)
    end
    filter(x->x.rank==1, P)
end

function fast_non_dominated_sort!(pop::AbstractVector{T}) where {T}
    F = T[]
    n = length(pop)

    for p in pop
        empty!(p.dom_list)
        p.dom_count = 0
    end

    for i in 1:n
        for j in i+1:n
            if pop[i] ⋖ pop[j]
                push!(pop[i].dom_list, j)
                pop[j].dom_count += 1
            elseif pop[j] ⋖ pop[i]
                push!(pop[j].dom_list, i)
                pop[i].dom_count += 1
            end
        end
        if pop[i].dom_count == 0
            pop[i].rank = 1
            push!(F, pop[i])
        end
    end
    res = Vector{T}[]
    i = 2
    while !isempty(F)
        Q = T[]
        for p in F
            for q in p.dom_list
                pop[q].dom_count -= 1
                if pop[q].dom_count == 0
                    pop[q].rank = i
                    push!(Q, pop[q])
                end
            end
        end
        i += 1
        push!(res, F)
        F = Q
    end
    res
end

function crowding_distance_assignement!(pop::AbstractVector{indiv{X,G,2,Y}}) where {X, G, Y}
    sort!(pop, by = x-> x.y[1])
    pop[1].crowding = pop[end].crowding = Inf
    for i = 2:length(pop)-1
        pop[i].crowding = (pop[i+1].y[1]-pop[i-1].y[1]) / (pop[end].y[1]-pop[1].y[1])
        pop[i].crowding += (pop[i-1].y[2]-pop[i+1].y[2]) / (pop[1].y[2]-pop[end].y[2])
    end
end

function crowding_distance_assignement!(pop::AbstractVector{indiv{X,G,N,Y}}) where {X,G,N,Y}
    for ind in pop
        ind.crowding = 0.
    end
    @inbounds for j = 1:N # Foreach objective
        let j = j #https://github.com/JuliaLang/julia/issues/15276
            sort!(pop, by = x-> x.y[j]) #sort by the objective value
        end
        pop[1].crowding = pop[end].crowding = Inf #Assign infinite value to extremas
        if pop[1].y[j] != pop[end].y[j]
            for i = 2:length(pop)-1
                pop[i].crowding += (pop[i+1].y[j]-pop[i-1].y[j]) / (pop[end].y[j]-pop[1].y[j])
            end
        end
    end
end

function tournament_selection(P)
    a, b = rand(P), rand(P)
    a < b ? a : b
end