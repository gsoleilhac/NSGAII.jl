# I.Y. Kim and O.L. de Weck.  Adaptive weighted-sum method for bi-objective optimization
# Pareto front generation. Structural and Multidisciplinary Optimization.
# Vol. 29, Num. 2, pp. 149–158, 2005.
# https://link.springer.com/article/10.1007/s00158-004-0465-1

using NSGAII
using Plots: scatter

function plot_pop(pop)
    pop = filter(x -> x.rank == 1, pop)
    display(scatter(map(x -> x.y[1], pop), map(x -> x.y[2], pop), markersize = 1))
    sleep(0.1)
end

const d = BinaryCoding(4, [-3, -3], [3, 3])
z1(x1, x2) = -(3(1-x1)^2 * exp(-x1^2 - (x2+1)^2) - 10(x1/5 - x1^3 - x2^5) * exp(-x1^2-x2^2) -3exp(-(x1+2)^2 - x2^2) + 0.5(2x1 + x2))
z2(x1, x2) = -(3(1+x2)^2 * exp(-x2^2 - (1-x1)^2) - 10(-x2/5 + x2^3 + x1^5) * exp(-x1^2-x2^2) - 3exp(-(2-x2)^2 - x1^2))
z(x) = z1(x[1], x[2]), z2(x[1], x[2])
nsga(300, 40, z, d, fplot = plot_pop, plotevery = 2)
