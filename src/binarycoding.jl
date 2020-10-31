struct BinaryCoding
    nbvar::Int
    types::Vector{Symbol}
    lb::Vector{Float64}
    ub::Vector{Float64}
    nbbits::Vector{Int}
    nbbitstotal::Int
end

struct BinaryCodedIndiv
    x::BitVector
    p::Vector{Float64}
end

function BinaryCoding(ϵ::Int, types::Vector{Symbol}, lb, ub)
    @assert length(types) == length(lb) == length(ub)
    @assert all(lb .< ub)
    @assert UInt128(10)^ϵ <= UInt128(2)^127

    nbvar = length(lb)
    nbbits = ones(Int, nbvar)
    for i = 1:nbvar
        if types[i] == :Cont
            while ((UInt128(10)^ϵ)*(ub[i]-lb[i]) >= UInt128(2)^nbbits[i])
                nbbits[i] += 1
            end
        elseif types[i] == :Int
            while 2^nbbits[i] <= ub[i]-lb[i]
                nbbits[i] += 1
            end
        end
    end
    BinaryCoding(nbvar, types, lb, ub, nbbits, sum(nbbits))
end
BinaryCoding(ϵ::Int, lb, ub) = BinaryCoding(ϵ, fill(:Cont, length(lb)), lb, ub)

function decode!(indiv::BinaryCodedIndiv, d::BinaryCoding)
    j = 0
    for i = 1:d.nbvar
        j += d.nbbits[i]
        if d.types[i] == :Bin
            indiv.p[i] = indiv.x[j] == 1 ? 1. : 0.
        else
            val = zero(UInt128)
            puis = one(UInt128)
            for ind = j:-1:j-d.nbbits[i]+1
                indiv.x[ind] && (val += puis)
                puis *= 2
            end

            if d.types[i] == :Cont
                indiv.p[i] = d.lb[i] + val * (d.ub[i] - d.lb[i]) / (UInt128(2)^d.nbbits[i] - 1)
            else
                indiv.p[i] = val + d.lb[i]
            end
        end
    end
    indiv
end

encode(p, d) = encode([p], d)
function encode(p::AbstractVector, d::BinaryCoding)::BinaryCodedIndiv
    res = BinaryCodedIndiv(BitVector(), p)
    sizehint!(res.x, d.nbbitstotal)
    for i = 1:d.nbvar
        if d.types[i] == :Int
            tab = reverse(digits(Bool, round(Int, res.p[i] - d.lb[i]), base = 2, pad = d.nbbits[i]))
            append!(res.x, tab)
        elseif d.types[i] == :Bin
            push!(res.x, p[i]!=0)
        else
            t = (p[i] - d.lb[i]) / (d.ub[i] - d.lb[i]) * (UInt128(2)^d.nbbits[i] - 1)
            target = round(UInt128, t)
            if target == UInt128(2)^d.nbbits[i] - 1
                target -= 1
            end
            tab = reverse(digits(Bool, target, base = 2, pad = d.nbbits[i]))
            append!(res.x, tab)
        end
    end
    res
end
