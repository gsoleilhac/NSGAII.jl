using NSGAII
using Base.Test

# write your own tests here
const C1 = [2 5 4 7 ; 3 3 5 7 ; 3 8 4 2 ; 6 5 2 5]
const C2 = [3 3 6 2 ; 5 3 7 3 ; 5 2 7 4 ; 4 6 3 5]
const C3 = [4 2 5 3 ; 5 3 4 3 ; 4 3 5 2 ; 6 4 7 3]

z(x, C) = sum(inds->C[inds...], enumerate(x))
z(x::Vector{Int}) = z(x, C1), z(x, C2), z(x, C3)

res = unique(nsga(500, 100, ()->randperm(4), z))
@show res
@test length(res) == 7


const d = RealData(6, [-10], [10])
z1(x) = x^2
z2(x) = (x-2)^2
z(bits) = begin 
    x = decode(bits, d)[1]
    z1(x), z2(x)
end
seed = encode.([-10 + rand()*20 for _ =1:100] ,d)
res = nsga(500, 200, ()->rand(Bool, d.nbbitstotal), z, seed = seed)

@test maximum(x -> x.y[1], res) >= 3.99
@test minimum(x -> x.y[1], res) <= 1e-4
@test maximum(x -> x.y[2], res) >= 3.99
@test minimum(x -> x.y[2], res) <= 1e-4
