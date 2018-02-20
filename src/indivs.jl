mutable struct indiv{G, P, N, Y<:Number}#Genotype, Phenotype, Dimension(Y_N), Type(Y_N)
    x::G
    pheno::P
    y::SVector{N, Y}
    CV::Float64
    rank::UInt16
    crowding::Float64
    dom_count::UInt16
    dom_list::Vector{UInt16}
    indiv(x::G, pheno::P, y::SVector{N, Y}, cv) where {G,P,N,Y} = new{G, P, N, Y}(x, pheno, y, cv, 0, 0., 0, [])
    indiv(x::G, pheno::P, y::NTuple{N, Y}, cv) where {G,P,N,Y} = indiv(x, pheno, SVector(y...), cv)
end


function indiv(geno, funcs::NSGA_FUNCS{G,P,N,Y}) where {G,P,N,Y}
    pheno = funcs.fdecode(convert(G, geno))
    indiv(geno, pheno, promote(funcs.feval(pheno)...), funcs.fCV(pheno))
end
indiv(funcs) = indiv(funcs.finit(), funcs)

struct Max end
struct Min end

function dominates(::Min, a::indiv, b::indiv)
    a.CV != b.CV && return a.CV < b.CV
    res = false
    for i in eachindex(a.y)
        @inbounds a.y[i] > b.y[i] && return false
        @inbounds a.y[i] < b.y[i] && (res=true)
    end
    res
end

function dominates(::Max, a::indiv, b::indiv)
    a.CV != b.CV && return a.CV < b.CV
    res = false
    for i in eachindex(a.y)
        @inbounds a.y[i] < b.y[i] && return false
        @inbounds a.y[i] > b.y[i] && (res=true)
    end
    res
end

Base.:(==)(a::indiv, b::indiv) = a.x == b.x
Base.hash(a::indiv) = hash(a.x)
Base.isless(a::indiv, b::indiv) = a.rank < b.rank || a.rank == b.rank && a.crowding >= b.crowding #Comparison operator for tournament selection
(Base.show(io::IO, ind::indiv{G,P,Y,N}) where {G,P,Y,N}) = print(io, "ind($(ind.pheno) : $(ind.y) | rank : $(ind.rank))")
(Base.show(io::IO, ind::indiv{G,P,N,Y}) where {G,P<:AbstractVector{Bool},N,Y}) = print(io, "ind($(String(map(x->x ? '1' : '0', ind.pheno))) : $(ind.y) | rank : $(ind.rank))")

function eval!(indiv::indiv, fdecode, z, fCV)
    indiv.pheno = fdecode(indiv.x)
    indiv.CV = fCV(indiv.pheno)
    indiv.CV ≈ 0 && (indiv.y = z(indiv.pheno))
    indiv
end