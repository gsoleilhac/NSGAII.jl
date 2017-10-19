function rand_flip!(bits)
    nb = length(bits)
    for i = 1:nb
        if rand() < 1/nb
            bits[i] = 1 - bits[i]
        end
    end
end

default_mutation!(b::BitArray) = rand_flip!(b)
default_mutation!(b::Vector{Bool}) = rand_flip!(b)

function rand_swap!(perm::Vector{Int})
    i = j = rand(1:length(perm))
    while j == i
        j = rand(1:length(perm))
    end
    perm[i], perm[j] = perm[j], perm[i]
end

default_mutation!(p::Vector{Int}) = rand_swap!(p)


mutate!(ind::indiv, f!) = f!(ind.x)