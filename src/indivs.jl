mutable struct indiv{X, Y}
    x::X
    y::Y
    rank::Int
    crowding::Float64
    dom_count::Int
    dom_list::Vector{indiv{X,Y}}
    indiv(x::X, y::Y) where {X,Y} = new{X, Y}(x, y, 0, 0., 0, indiv{X,Y}[])
end
indiv(x, z::Function) = indiv(x, z(x))


function dominates(a::indiv, b::indiv)
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
import Base.show; show(io::IO, ind::indiv) = print(io, "ind($(ind.x) : $(ind.y) | rank : $(ind.rank))")
function show(io::IO, ind::indiv{X, Y}) where {X<:AbstractArray{Bool}, Y}
    print(io, "ind(") ; show(IOContext(io, limit=true, compact=true), Int.(ind.x)) ; print(io, " : $(ind.y) | rank : $(ind.rank))")
end

function eval!(indiv::indiv, z::Function)
    indiv.y = z(indiv.x)
    indiv
end