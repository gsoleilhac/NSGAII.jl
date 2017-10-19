# G. Mavrotas, D. Diakoulaki.  
# A branch and bound algorithm for mixed zero-one multiple objective linear programming.
# European Journal of Operational Research, 107, pp. 530–541, 1998.
# http://www.sciencedirect.com/science/article/pii/S0377221797000775


using NSGAII, PyPlot

function plot_pop(P)
    clf()
    pop = filter(x -> all(x.y .< 5000), P)
    non_dom = NSGAII.fast_non_dominated_sort!(P)[1]
    pop = setdiff(pop, non_dom)
    p = plot3D(map(x -> x.y[1], pop), map(x -> x.y[2], pop), map(x -> x.y[3], pop), "bo", markersize=1)
    p = plot3D(map(x -> x.y[1], non_dom), map(x -> x.y[2], non_dom), map(x -> x.y[3], non_dom), "go", markersize=1)
    ax = gca()
    ax[:set_xlim]([-95, 110])
    ax[:set_ylim]([-60, 45])
    ax[:set_zlim]([-50, 110])
    !isinteractive() && show()
    sleep(0.1)
end

z1(x, δ) = dot([17,-12,-12,-19,-6], x) + dot([-73, -99, -81], δ)
z2(x, δ) = dot([2,-6,0,-12,13], x) + dot([-61,-79,-53], δ)
z3(x, δ) = dot([-20,7,-16,0,-1], x) + dot([-72,-54,-79], δ)


d = RealCoding(6, zeros(5), fill(10, 5))
init = ()->rand(Bool, d.nbbitstotal + 3) #+ 3 bits pour δ
_decode(x) = decode(x[1:end-3], d), x[end-2:end] #renvoie un tuple : Tableau de 5 float décodés, + tableau des 3 derniers booléens tels quels.

function isvalid(x, δ)
    sum(δ) > 1 && return false
    -x[2] + 6x[5] + 25δ[1] > 52 && return false
    -x[1] + 18x[4] + 18x[5] + 8δ[2] > 77 && return false
    7x[4] + 9x[5] + 19δ[3] > 66 && return false
    16x[1] + 20x[5] > 86 && return false
    13x[2] + 7x[4] > 86 && return false
    true
end

function z_ex4(bits)
    x, δ = _decode(bits)
    !isvalid(x, δ) && return (5000.,5000.,5000.)
    return -z1(x, δ), -z2(x, δ), -z3(x, δ)
end

lex_1 = vcat( encode([5.375, 0., 0., 0., 0.], d), false, false, false )
lex_2 = vcat( encode([0.025974, 0., 0., 0., 4.27922], d), false, false, false )
lex_3 = vcat( encode([0., 6.61538, 0., 0., 0.], d), false, false, false )

@time pop = nsga(500, 50, init, z_ex4, seed=[lex_1, lex_2, lex_3] , fplot = plot_pop )
@time pop = nsga(200, 500, init, z_ex4)
for i in pop
    i.y = (-).(i.y)
end
