function fast_non_dominated_sort!(pop::Vector{T}) where {T}
    F = T[]
    for p in pop
        empty!(p.dom_list)
        p.dom_count = 0
        for q in pop
            if p ⋖ q
                push!(p.dom_list, q)
            elseif q ⋖ p
                p.dom_count += 1
            end
        end
        if p.dom_count == 0
            p.rank = 1
            push!(F, p)
        end
    end
    res = Vector{T}[]
    i = 2
    while !isempty(F)
        Q = T[]
        for p in F
            for q in p.dom_list
                q.dom_count -= 1
                if q.dom_count == 0
                    q.rank = i
                    push!(Q, q)
                end
            end
        end
        i += 1
        push!(res, F)
        F = Q
    end
    res
end

function crowding_distance_assignement!(pop::Vector{indiv{X,2,Y}}) where {X, Y}
    sort!(pop, by = x-> x.y[1])
    pop[1].crowding = pop[end].crowding = Inf
    for i = 2:length(pop)-1
        pop[i].crowding = (pop[i+1].y[1]-pop[i-1].y[1]) / (pop[end].y[1]-pop[1].y[1])
        pop[i].crowding += (pop[i-1].y[2]-pop[i+1].y[2]) / (pop[1].y[2]-pop[end].y[2])
    end
end

function crowding_distance_assignement!(pop::Vector{indiv{X,N,Y}}) where {X,N,Y}
    for ind in pop
        ind.crowding = 0.
    end
    for j = 1:length(pop[1].y) # Foreach objective
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