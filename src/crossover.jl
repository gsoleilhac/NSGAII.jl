function crossover!(ind_a, ind_b, fcross, child_a, child_b)
    fcross(ind_a.x, ind_b.x, child_a.x, child_b.x)
end
(default_crossover!(pa::T, pb::T, ca, cb)) where T<:AbstractVector{Bool} = two_point_crossover!(pa, pb, ca, cb)
(default_crossover!(pa::T, pb::T, ca, cb)) where T<:AbstractVector{Int} = PMX_crossover!(pa, pb, ca, cb)

function two_point_crossover!(bits_a, bits_b, child1, child2)
    cut_a = cut_b = rand(2:length(bits_a)-1)
    while(cut_b == cut_a)
        cut_b = rand(2:length(bits_a))
    end
    cut_a,cut_b = minmax(cut_b,cut_a)

    copyto!(child1, 1, bits_a, 1, cut_a-1)
    copyto!(child1, cut_a, bits_b, cut_a, cut_b-cut_a+1)
    copyto!(child1, cut_b+1, bits_a, cut_b+1, length(bits_a)-cut_b)

    copyto!(child2, 1, bits_b, 1, cut_a-1)
    copyto!(child2, cut_a, bits_a, cut_a, cut_b-cut_a+1)
    copyto!(child2, cut_b+1, bits_b, cut_b+1, length(bits_a)-cut_b)
end

function PMX_crossover!(pa, pb, ca, cb)
    cut_a = cut_b = rand(1:length(pa))
    while(cut_b == cut_a)
        cut_b = rand(1:length(pa))
    end
    cut_a, cut_b = minmax(cut_a, cut_b)

    copyto!(ca, 1, pb, 1, cut_a-1)
    copyto!(ca, cut_a, pa, cut_a, cut_b-cut_a+1)
    copyto!(ca, cut_b+1, pb, cut_b+1, length(pa)-cut_b)

    copyto!(cb, 1, pa, 1, cut_a-1)
    copyto!(cb, cut_a, pb, cut_a, cut_b-cut_a+1)
    copyto!(cb, cut_b+1, pa, cut_b+1, length(pa)-cut_b)

    @inbounds for i = cut_a:cut_b
        if pa[i] ∉ view(pb, cut_a:cut_b)
            j = findfirst(isequal(pb[i]), pa)
            while j ∈ cut_a:cut_b
                j = findfirst(isequal(pb[j]), pa)
            end
            cb[j] = pa[i]
        end

        if pb[i] ∉ view(pa, cut_a:cut_b)
            j = findfirst(isequal(pa[i]), pb)
            while j ∈ cut_a:cut_b
                j = findfirst(isequal(pa[j]), pb)
            end
            ca[j] = pb[i]
        end
    end
end