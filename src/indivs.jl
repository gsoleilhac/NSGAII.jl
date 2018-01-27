mutable struct indiv{G, P, N, Y<:Number}#Genotype, Phenotype, Dimension(Y_N), Type(Y_N)
    x::G
    pheno::P
    y::NTuple{N, Y}
    CV::Float64
    rank::Int
    crowding::Float64
    dom_count::Int
    dom_list::Vector{indiv{G,P,N,Y}}
    indiv(x::G, pheno::P, y::NTuple{N, Y}, cv) where {G,P,N,Y} = new{G, P, N, Y}(x, pheno, y, cv, 0, 0., 0, indiv{G,P,N,Y}[])
end
function indiv(x, fdecode::Function, z::Function, fCV::Function)
    pheno = fdecode(x)
    indiv(x, pheno, z(pheno), fCV(pheno))
end

function dominates(a::indiv, b::indiv)
    a.CV != b.CV && return a.CV < b.CV

    res = false
    for i in eachindex(a.y)
        a.y[i] > b.y[i] && return false
        a.y[i] < b.y[i] && (res=true)
    end
    res
end
⋖(a::indiv, b::indiv) = dominates(a, b)

import Base.==; ==(a::indiv, b::indiv) = a.x == b.x
import Base.hash; hash(a::indiv) = hash(a.x)
import Base.isless; isless(a::indiv, b::indiv) = a.rank < b.rank || a.rank == b.rank && a.crowding >= b.crowding #Comparison operator for tournament selection
import Base.show; show(io::IO, ind::indiv) = print(io, "ind($(ind.pheno) : $(ind.y) | rank : $(ind.rank))")

function eval!(indiv::indiv, fdecode::Function, z::Function, fCV::Function)
    indiv.pheno = fdecode(indiv.x)
    indiv.CV = fCV(indiv.pheno)
    indiv.CV ≈ 0 && (indiv.y = z(indiv.pheno))
    indiv
end