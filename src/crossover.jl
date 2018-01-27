function crossover(ind_a, ind_b, fcross)
    c1, c2 = fcross(ind_a.x, ind_b.x)
    return indiv(c1, ind_a.pheno, ind_a.y, 0.), indiv(c2, ind_a.pheno, ind_a.y, 0.) #phenotypes, objective value and CV value will be calculated later, after mutation
end

function two_point_crossover(bits_a, bits_b)
    cut_a = cut_b = rand(2:length(bits_a)-1)
    while(cut_b == cut_a)
        cut_b = rand(2:length(bits_a))
    end
    cut_a,cut_b = minmax(cut_b,cut_a)
    child1 = vcat(bits_a[1:cut_a-1], bits_b[cut_a:cut_b], bits_a[cut_b+1:end])
    child2 = vcat(bits_b[1:cut_a-1], bits_a[cut_a:cut_b], bits_b[cut_b+1:end])
    child1, child2
end
(default_crossover(ba::T, bb::T)) where T<:AbstractVector{Bool} = two_point_crossover(ba, bb)

function PMX_crossover(pa, pb)
    cut_a = cut_b = rand(1:length(pa))
    while(cut_b == cut_a)
        cut_b = rand(1:length(pa))
    end
    cut_a, cut_b = minmax(cut_a, cut_b)
    ca, cb = copy(pb), copy(pa)
    ca[cut_a:cut_b] .= pa[cut_a:cut_b]
    cb[cut_a:cut_b] .= pb[cut_a:cut_b]

    for i = cut_a:cut_b
        if pa[i] ∉ pb[cut_a:cut_b]
            j = findfirst(pa, pb[i])
            while j ∈ cut_a:cut_b
                j = findfirst(pa, pb[j])
            end
            cb[j] = pa[i]
        end

        if pb[i] ∉ pa[cut_a:cut_b]
            j = findfirst(pb, pa[i])
            while j ∈ cut_a:cut_b
                j = findfirst(pb, pa[j])
            end
            ca[j] = pb[i]
        end
    end
    ca, cb
end
(default_crossover(ba::T, bb::T)) where T<:AbstractVector{Int} = PMX_crossover(ba, bb)