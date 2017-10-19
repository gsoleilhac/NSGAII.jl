# Anthony Przybylski, Xavier Gandibleux, Matthias Ehrgott.
# Two phase algorithms for the bi-objective assignment problem.
# European Journal of Operational Research, Vol. 185, Iss. 2, pp. 509–533, 2008.
# http://www.sciencedirect.com/science/article/pii/S0377221707000781

# Anthony Przybylski, Xavier Gandibleux, Matthias Ehrgott.
# A two phase method for multi-objective integer programming 
# and its application to the assignment problem with three objectives.
# Discrete Optimization, Vol. 7, Iss. 3, pp. 149–165, 2010.
# http://www.sciencedirect.com/science/article/pii/S1572528610000125

using NSGAII

const C1 = [2 5 4 7 ; 3 3 5 7 ; 3 8 4 2 ; 6 5 2 5]
const C2 = [3 3 6 2 ; 5 3 7 3 ; 5 2 7 4 ; 4 6 3 5]
const C3 = [4 2 5 3 ; 5 3 4 3 ; 4 3 5 2 ; 6 4 7 3]

z(x, C) = sum(inds->C[inds...], enumerate(x))
z(x::Vector{Int}) = z(x, C1), z(x, C2), z(x, C3)

res = nsga(200, 50, ()->randperm(4), z)

for x in unique(res)
    println(x)
end
