struct MixedCoding
    nbvar::Int
    types::Vector{Symbol}
    lb::Vector{Float64}
    ub::Vector{Float64}
    nbbits::Vector{Int}
    nbbitstotal::Int
end

function MixedCoding(ϵ::Int, types::Vector{Symbol}, lb, ub)
    @assert length(types)==length(lb)==length(ub)
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
    MixedCoding(nbvar, types, lb, ub, nbbits, sum(nbbits))
end
RealCoding(ϵ::Int, lb, ub) = MixedCoding(ϵ, fill(:Cont, length(lb)), lb, ub)

function decode(x, d::MixedCoding)::Vector{Float64}

    res = zeros(d.nbvar)
    for i = 1:d.nbvar
        if d.types[i] == :Bin
            res[i] = x[sum(d.nbbits[1:i])]==1 ? 1. : 0.
        else
            val = zero(UInt128)
            puis = one(UInt128)
            jstart = sum(d.nbbits[1:i])
            jend = jstart - d.nbbits[i] + 1
            for j = jstart:-1:jend
                x[j] && (val += puis)
                puis *= 2
            end

            if d.types[i] == :Cont
                res[i] = d.lb[i] + val * (d.ub[i] - d.lb[i]) / (UInt128(2)^d.nbbits[i] - 1)
            else
                res[i] = val - d.lb[i]
            end
        end
    end
    res
end

function encode(x, d::MixedCoding)::BitVector
    res = BitVector(0)
    sizehint!(res, d.nbbitstotal)
    for i = 1:d.nbvar
        if d.types[i] == :Int
            tab = reverse(digits(Bool, round(Int, x[i] + d.lb[i]), 2, d.nbbits[i]))
            append!(res, tab)
        elseif d.types[i] == :Bin
            push!(res, x[i]!=0)
        else
            t = (x[i] - d.lb[i]) / (d.ub[i] - d.lb[i]) * (UInt128(2)^d.nbbits[i] - 1)
            @show t
            target = round(UInt128, t)
            if target == UInt128(2)^d.nbbits[i] - 1
                target -= 1
            end
            tab = reverse(digits(Bool, target, 2, d.nbbits[i]))
            append!(res, tab)
        end
    end
    res
end
