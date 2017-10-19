struct RealData
    nbvar::Int
    lb::Vector{Float64}
    ub::Vector{Float64}
    nbbits::Vector{Int}
    nbbitstotal::Int
end
RealData(ϵ::Vector{Int}, lb, ub) = begin
    nbvar = length(lb)
    nbbits = ones(Int, nbvar)
    for i = 1:nbvar
        while ((ub[i]-lb[i])*(10^ϵ[i]) >= 2^nbbits[i])
            nbbits[i] += 1
        end
    end
    RealData(nbvar, lb, ub, nbbits, sum(nbbits))
end
RealData(ϵ::Int, lb, ub) = RealData([ϵ for _=1:length(lb)], lb, ub)

function decode(x, d::RealData)::Vector{Float64}

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

function encode(x, d::RealData)::Vector{Bool}

    res = Bool[]
    for i = 1:d.nbvar
        tab = [false for i = 1:d.nbbits[i]]
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