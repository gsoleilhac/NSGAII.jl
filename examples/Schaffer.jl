# David Schaffer. Multiple Objective Optimization with Vector Evaluated Genetic Algorithms.
# In Proceedings of the 1st International Conference on Genetic Algorithms, 
# L. Erlbaum Associates Inc. pp. 93–100, 1985.
# http://dl.acm.org/citation.cfm?id = 645511.657079

using NSGAII, Plots

function plot_pop(pop)
    pop = filter(x -> x.rank == 1, pop)
    display(scatter(map(x -> x.y[1], pop), map(x -> x.y[2], pop), markersize = 1))
    sleep(0.1)
end

const d = BinaryCoding(6, [-10], [10])
z(x) = ( x[1]^2 , (x[1]-2)^2 )
nsga(200, 30, z, d, fplot = plot_pop)
