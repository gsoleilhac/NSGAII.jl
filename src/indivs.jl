mutable struct Indiv{G, Y} # Genotype, Type(Y_N)
    x::G
    y::Y
    CV::Float64
    rank::UInt16
    crowding::Float64
    dom_count::UInt16
    dom_list::Vector{UInt16}
    Indiv(x::G, y::Y, cv) where {G, Y} = new{G, Y}(x, y, cv, zero(UInt16), 0., zero(UInt16), UInt16[])
end
function createIndiv(x, z, fCV)
    y = z(x)
    cv = fCV(x)
    Indiv(x, y, cv)
end

struct Max end
struct Min end

function dominates(::Min, a::Indiv, b::Indiv)
    a.CV != b.CV && return a.CV < b.CV
    res = false
    for i in eachindex(a.y)
        @inbounds a.y[i] > b.y[i] && return false
        @inbounds a.y[i] < b.y[i] && (res = true)
    end
    res
end

function dominates(::Max, a::Indiv, b::Indiv)
    a.CV != b.CV && return a.CV < b.CV
    res = false
    for i in eachindex(a.y)
        @inbounds a.y[i] < b.y[i] && return false
        @inbounds a.y[i] > b.y[i] && (res = true)
    end
    res
end

Base.:(==)(a::Indiv, b::Indiv) = a.x == b.x
Base.hash(a::Indiv) = hash(a.x)
Base.isless(a::Indiv, b::Indiv) = a.rank < b.rank || a.rank == b.rank && a.crowding >= b.crowding #Comparison operator for tournament selection
Base.show(io::IO, ind::Indiv) = print(io, "Indiv($(repr_pheno(ind.x)) : $(ind.y) | rank : $(ind.rank))")
repr_pheno(x) = repr(x)
function repr_pheno(x::Union{BitVector, Vector{Bool}}) 
    res = map(x -> x ? '1' : '0', x)
    if length(res) <= 40
        return "["*String(res)*"]"
    else
        return "["*String(res[1:15])*"..."*String(res[end-14:end])*"]"
    end
end

function eval!(indiv::Indiv, z::Function, fCV::Function)
    indiv.CV = fCV(indiv.x)
    indiv.CV ≈ 0 && (indiv.y = z(indiv.x))
    indiv
end