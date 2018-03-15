# I.Y. Kim and O.L. de Weck.  Adaptive weighted-sum method for bi-objective optimization
# Pareto front generation. Structural and Multidisciplinary Optimization.
# Vol. 29, Num. 2, pp. 149–158, 2005.
# https://link.springer.com/article/10.1007/s00158-004-0465-1

using NSGAII, PyPlot

function plot_pop(P)
    clf()
    p = plot(map(x -> x.y[1], P), map(x -> x.y[2], P), "bo", markersize=1)
    !isinteractive() && show()
    sleep(0.2)
end

const d = BinaryCoding(8, [-3, -3], [3, 3])
z1(x1, x2) = -(3(1-x1)^2 * exp(-x1^2 - (x2+1)^2) - 10(x1/5 - x1^3 - x2^5) * exp(-x1^2-x2^2) -3exp(-(x1+2)^2 - x2^2) + 0.5(2x1 + x2))
z2(x1, x2) = -(3(1+x2)^2 * exp(-x2^2 - (1-x1)^2) - 10(-x2/5 + x2^3 + x1^5) * exp(-x1^2-x2^2) - 3exp(-(2-x2)^2 - x1^2))
z(x) = z1(x[1], x[2]), z2(x[1], x[2])
nsga(300, 20, z, d, fplot = plot_pop)
