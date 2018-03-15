using NSGAII, vOptGeneric, vOptSpecific
using Compat.Test, Compat.Random



const C1 = [2 5 4 7 ; 3 3 5 7 ; 3 8 4 2 ; 6 5 2 5]
const C2 = [3 3 6 2 ; 5 3 7 3 ; 5 2 7 4 ; 4 6 3 5]
const C3 = [4 2 5 3 ; 5 3 4 3 ; 4 3 5 2 ; 6 4 7 3]
z(x, C) = sum(inds->C[inds...], enumerate(x))
z(x::Vector{Int}) = z(x, C1), z(x, C2), z(x, C3)
res = unique(nsga(500, 100, z, ()->randperm(4)))
@test length(res) == 7

const d = RealCoding(6, [-10], [10])
z(x) = x[1]^2, (x[1] - 2)^2
seed = [-10 + rand()*20 for _ =1:100]
res = nsga(500, 200, z, d, seed = seed)
@test maximum(x -> x.y[1], res) >= 3.99
@test minimum(x -> x.y[1], res) <= 1e-4
@test maximum(x -> x.y[2], res) >= 3.99
@test minimum(x -> x.y[2], res) <= 1e-4

m = vModel()
id = load2UKP(joinpath(@__DIR__, "2KP500-1A.DAT"))
p1, p2, w, c = id.P1, id.P2, id.W, id.C
@variable(m, x[1:length(p1)], Bin)
@addobjective(m, Max, dot(x, p1))
@addobjective(m, Max, dot(x, p2))
@constraint(m, dot(x, w) <= c)
@inferred nsga_binary(100, 200, m)

m = vModel()
@variable(m, 0 <=x[1:5] <= 10)
@variable(m, 0 <= δ[1:3] <= 1, Int)
@addobjective(m, Max, dot([17,-12,-12,-19,-6], x) + dot([-73, -99, -81], δ))
@addobjective(m, Max, dot([2,-6,0,-12,13], x) + dot([-61,-79,-53], δ))
@addobjective(m, Max, dot([-20,7,-16,0,-1], x) + dot([-72,-54,-79], δ))
@constraint(m, sum(δ) <= 1)
@constraint(m, -x[2] + 6x[5] + 25δ[1] <= 52)
@constraint(m, -x[1] + 18x[4] + 18x[5] + 8δ[2] <= 77)
@constraint(m, 7x[4] + 9x[5] + 19δ[3] <= 66)
@constraint(m, 16x[1] + 20x[5] <= 86)
@constraint(m, 13x[2] + 7x[4] <= 86)
@inferred nsga(100, 200, m)