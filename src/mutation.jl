mutate!(ind::Indiv, fmut!) = fmut!(ind.x)
default_mutation!(p::Vector{Int}) = rand_swap!(p)
default_mutation!(b::T) where T<:AbstractVector{Bool} = rand_flip!(b)

function rand_flip!(bits)
    nb = length(bits)
    for i = 1:nb
        if rand() < 1/nb
            @inbounds bits[i] = 1 - bits[i]
        end
    end
end

function rand_swap!(perm::Vector{Int})
    i = j = rand(1:length(perm))
    while j == i
        j = rand(1:length(perm))
    end
    @inbounds perm[i], perm[j] = perm[j], perm[i]
end

