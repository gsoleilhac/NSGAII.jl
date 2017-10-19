# David Schaffer. Multiple Objective Optimization with Vector Evaluated Genetic Algorithms.
# In Proceedings of the 1st International Conference on Genetic Algorithms, 
# L. Erlbaum Associates Inc. pp. 93–100, 1985.
# http://dl.acm.org/citation.cfm?id=645511.657079

using NSGAII, PyPlot

function plot_pop(P)
    clf()
    p = plot(map(x -> x.y[1], P), map(x -> x.y[2], P), "bo", markersize=1)
    !isinteractive() && show()
    sleep(0.2)
end

const d = RealCoding(6, [-10], [10])
z1(x) = x^2
z2(x) = (x-2)^2
z(bits) = begin 
    x = decode(bits, d)[1]
    z1(x), z2(x)
end
nsga(200, 20, ()->rand(Bool, d.nbbitstotal), z, plot_func = plot_pop)
