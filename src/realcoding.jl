struct RealCoding
    nbvar::Int
    lb::Vector{Float64}
    ub::Vector{Float64}
    nbbits::Vector{Int}
    nbbitstotal::Int
end
RealCoding(ϵ::Vector{Int}, lb, ub) = begin
    @assert length(lb)==length(ub)
    for i = 1:length(lb)
        @assert lb[i] < ub[i]
    end
    nbvar = length(lb)
    nbbits = ones(Int, nbvar)
    for i = 1:nbvar
        while ((Int128(10)^ϵ[i])*(ub[i]-lb[i]) >= Int128(2)^nbbits[i])
            nbbits[i] += 1
        end
    end
    RealCoding(nbvar, lb, ub, nbbits, sum(nbbits))
end
RealCoding(ϵ::Int, lb, ub) = RealCoding([ϵ for _=1:length(lb)], lb, ub)

function decode(x, d::RealCoding)::Vector{Float64}

    res = zeros(d.nbvar)
    for i = 1:d.nbvar
        val = 0
        puis = 1
        jstart = sum(d.nbbits[1:i])
        jend = jstart - d.nbbits[i] + 1
        for j = jstart:-1:jend
            x[j] && (val += puis)
            puis *= 2
        end
        res[i] = d.lb[i] + val * (d.ub[i] - d.lb[i]) / (2^d.nbbits[i] - 1)
    end
    res
end

function encode(x, d::RealCoding)::BitVector

    res = BitArray(0)
    for i = 1:d.nbvar
        tab = falses(d.nbbits[i])
        ind = 1
        puis = 2^(d.nbbits[i] - 1)
        val = 0
        target = round(Int, (x[i] - d.lb[i]) * (2^d.nbbits[i] - 1) / (d.ub[i] - d.lb[i]))
        while val < target
            if val + puis <= target
                val += puis
                tab[ind] = true
            end
            puis ÷= 2
            ind += 1
        end
        append!(res, tab)
    end
    res
end