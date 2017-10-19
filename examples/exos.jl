using NSGAII
using PyPlot


function pyplot_ex4(P)
    clf()
    pop = filter(x -> all(x.y .< 5000), P)
    non_dom = fast_non_dominated_sort!(P)[1]
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


#EX 1
function ex1()
    d = RealData(6, [-10], [10])
    z_ex1(bits) = begin 
        x = decode(bits, d)[1]
        x^2, (x-2)^2
    end
    nsga(200, 20, ()->rand(Bool, d.nbbitstotal), z_ex1)
end

#EX 2
function ex2()
    z1_ex2(x1, x2) = -(3(1-x1)^2 * exp(-x1^2 - (x2+1)^2) - 10(x1/5 - x1^3 - x2^5) * exp(-x1^2-x2^2) -3exp(-(x1+2)^2 - x2^2) + 0.5(2x1 + x2))
    z2_ex2(x1, x2) = -(3(1+x2)^2 * exp(-x2^2 - (1-x1)^2) - 10(-x2/5 + x2^3 + x1^5) * exp(-x1^2-x2^2) - 3exp(-(2-x2)^2 - x1^2))
    d = RealData(8, [-3, -3], [6, 6])
    z_ex2(x::Vector{Bool}) = begin 
        x1, x2 = decode(x, d)
        z1_ex2(x1, x2), z2_ex2(x1, x2)
    end
    nsga(300, 20, ()->rand(Bool, d.nbbitstotal), z_ex2, plot_func = plot_func)
end

#EX 3
function ex3()
    C1 = [2 5 4 7 ; 3 3 5 7 ; 3 8 4 2 ; 6 5 2 5]
    C2 = [3 3 6 2 ; 5 3 7 3 ; 5 2 7 4 ; 4 6 3 5]
    C3 = [4 2 5 3 ; 5 3 4 3 ; 4 3 5 2 ; 6 4 7 3]
    z_ex3(x, C) = sum(inds->C[inds...], enumerate(x))
    z_ex3(x::Vector{Int}) = z_ex3(x, C1), z_ex3(x, C2), z_ex3(x, C3)
    nsga(200, 20, ()->randperm(4), z_ex3)
end

#EX 4
function ex4()
    z1(x, δ) = dot([17,-12,-12,-19,-6], x) + dot([-73, -99, -81], δ)
    z2(x, δ) = dot([2,-6,0,-12,13], x) + dot([-61,-79,-53], δ)
    z3(x, δ) = dot([-20,7,-16,0,-1], x) + dot([-72,-54,-79], δ)

    
    d = RealData(6, zeros(5), fill(10, 5))
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

    @time pop = nsga(500, 50, init, z_ex4, seed=[lex_1, lex_2, lex_3] ,plot_func = pyplot_ex4 )
    @time pop = nsga(200, 500, init, z_ex4)
    for i in pop
        i.y = (-).(i.y)
    end
    sort!(pop, by = ind -> ind.y[1])
    @show pop[end].y
    @show pop[1].y
    sort!(pop, by = ind -> ind.y[2])
    @show pop[end].y
    @show pop[1].y
    sort!(pop, by = ind -> ind.y[3])
    @show pop[end].y
    @show pop[1].y
end

function test()
    p = ex4()
    @profile ex4()
    ProfileView.view()
end

function test2()
    z1_ex2(x1, x2) = -(3(1-x1)^2 * exp(-x1^2 - (x2+1)^2) - 10(x1/5 - x1^3 - x2^5) * exp(-x1^2-x2^2) -3exp(-(x1+2)^2 - x2^2) + 0.5(2x1 + x2))
    z2_ex2(x1, x2) = -(3(1+x2)^2 * exp(-x2^2 - (1-x1)^2) - 10(-x2/5 + x2^3 + x1^5) * exp(-x1^2-x2^2) - 3exp(-(2-x2)^2 - x1^2))
    d = RealData(8, [-3, -3], [6, 6])
    z_ex2(x::Vector{Bool}) = begin 
        x1, x2 = decode(x, d)
        z1_ex2(x1, x2), z2_ex2(x1, x2)
    end
    P = [indiv(rand(Bool, d.nbbitstotal), z_ex2) for i=1:500]
    sort!(P, lt = (x,y) -> x.y[2] < y.y[2])
end