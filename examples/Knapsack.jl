using NSGAII, vOptGeneric, PyPlot
using GLPK, GLPKMathProgInterface ; solver = GLPKSolverMIP
# using CPLEX ; solver = CplexSolver()

n = 100
p1 = rand(60:100, n)
p2 = rand(60:100, n)
w = rand(50:100, n)
c = sum(w) ÷ 2

m = vModel(solver = GLPKSolverMIP())
@variable(m, x[1:length(p1)], Bin)
@addobjective(m, Max, dot(x, p1))
@addobjective(m, Max, dot(x, p2))
@constraint(m, dot(x, w) <= c)

function plot_pop(P)
    clf()
    P = filter(x -> x.rank == 1, P)
    plot(map(x -> x.y[1], P), map(x -> x.y[2], P), "bo", markersize = 2, label="nsga")
    show()
    sleep(0.1)
end

function greedy(p, w, c)
    order = sortperm(p./w, rev = true)
    sum_weight = 0
    res = falses(length(w))
    for i = 1:length(p)
        if w[order[i]] + sum_weight <= c
            res[order[i]] = true
            sum_weight += w[order[i]]
        end
    end
    res
end

seed = [greedy(α.*p1 .+ (1-α)*p2, w, c) for α=0:0.2:1]
nsga(100, 10000, m, fplot = plot_pop, seed = seed, plotevery = 1000)

solve(m, method=:epsilon)
Y_N = getY_N(m)
for n = 1:length(Y_N)
    X = getvalue(x, n)
    print(find(X))
    println("| z = ",Y_N[n])
end

f1 = map(y -> y[1], Y_N)
f2 = map(y -> y[2], Y_N)
xlabel("z1") ; ylabel("z2")
plot(f1,f2,"kx", markersize = 2, label="exact(ϵ-method)")
legend() ; show()

