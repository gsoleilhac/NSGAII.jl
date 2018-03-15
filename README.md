# Installation

```
Pkg.clone("https://github.com/gsoleilhac/NSGAII.jl")
```

# Usage

## Example : Bi-Objective Knapsack

```julia
using NSGAII
p1 = [77,94,71,63,96,82,85,75,72,91,99,63,84,87,79,94,90,60,69,62]
p2 = [65,90,90,77,95,84,70,94,66,92,74,97,60,60,65,97,93,60,69,74]
w = [80,87,68,72,66,77,99,85,70,93,98,72,100,89,67,86,91,79,71,99]
c = 900
n = length(p1)
```

The four mandatory parameters of NSGAII are 
* the size of the population
* the number of generations
* an evaluation function
* an initialization function

```julia
popsize = 100
nbgen = 200
init() = bitrand(n) #our genotype is a random binary vector
z(x) = dot(x, p1), dot(x, p2) #and our objectives are the sum of the items we pick
```
Now, this would be enough to run nsgaii with
`nsga_max(popsize, nbgen, z, init)`  
But we need to add the constraint that all items must fit in the knapsack.  
For this we define a *constraint-violation function* that returns 0 only if the solution is feasible,  
and a value > 0 otherwise.

```julia
function CV(x)
    sumW = dot(x, w)
    if sumW < c
        return 0
    else
        return sumW - c
    end
end

#We can now call
result = nsga_max(popsize, nbgen, z, init, fCV = CV)
```

### Crossover

Because the solutions are encoded as bitstrings, nsga will use by default a 2-point crossover, but we can define our own and assign it with the keyword `fcross`:

```julia
function one_point_crossover!(parent_a, parent_b, child_a, child_b)
    n = length(parent_a)
    cut = rand(1:n-1)

    child_a[1:cut] .= parent_a[1:cut]
    child_a[cut+1:n] .= parent_b[cut+1:n]

    child_b[1:cut] .= parent_b[1:cut]
    child_b[cut+1:n] .= parent_a[cut+1:n]
end

nsga_max(popsize, nbgen, z, init, fCV = CV, fcross = one_point_crossover!)
```
### Mutation

The default mutation for a binary vector is the *bitstring mutation* where each bit has a probability 1/l to be flipped (where l is the length of the vector)

As with crossovers, we can define or own mutation operator and assign it with the keyword `fmut`. The probability of mutation can be changed with the keyword `pmut`.

Let's say we want our mutation to flip two random bits :

```julia
function two_bits_flip!(bits)
    for i = 1:2
        n = rand(1:length(bits))
        bits[n] = 1 - bits[n]
    end
end

nsga_max(popsize, nbgen, z, init, fCV = CV, fmut = two_bits_flip!, pmut = 0.2)
```

### Genotype and Phenotype

So far, we haven't seen the difference between the genotype and the phenotype ; the default decoding function used here is the *identity* function

You can provide your own with the keywords `fdecode` and `fdecode!` which will work in-place.

Note : if your decode function takes a genotype `G` and returns a phenotype `P`, make sure your crossovers and mutations functions work on type `G`, and that your evaluation and (if provided) your constraint-violation  functions work on type `P`.  
`fdecode!` should take as parameters a genotype `G` and a phenotype `P` and modify it in-place.

See [BinaryCoding](https://github.com/gsoleilhac/NSGAII.jl#binarycoding)


### Seeding

Starting solutions can be provided as a vector with the keyword `seed`, for example : 

```julia
x1 = greedy(p1, w, c)
x2 = greedy(p2, w, c)
x3 = greedy(p1 .+ p2, w, c)

nsga_max(popsize, nbgen, z, init, fCV = CV, seed = [x1, x2, x3])
```

Make sure the type of your seeds is the same as the one given by calling `init()` !

### Plotting

A plot function can be passed with the keyword `fplot`, by default the population is plotted at every generation but this can be changed with the keyword `plotevery`.

Example with PyPlot : 
```julia
using PyPlot

function plot_pop(P)
    clf() #clears the figure
    P = filter(indiv -> indiv.rank == 1, P) #keep only the non-dominated solutions
    plot(map(x -> x.y[1], P), map(x -> x.y[2], P), "bo", markersize=1)
    sleep(0.1)
end

nsga_max(popsize, nbgen, z, init, fCV = CV, fplot = plot_pop, plotevery = 5)
```

## BinaryCoding

You can use `BinaryCoding(ϵ::Int, lb::Vector, ub::Vector)` to encode real variables with a precision `ϵ`, and with lower and upper bounds `lb` and `ub`

Example : 

```julia 
using NSGAII, PyPlot

function plot_pop(P)
    clf()
    P = filter(indiv -> indiv.rank == 1, P) #keep only the non-dominated solutions
    plot(map(x -> x.y[1], P), map(x -> x.y[2], P), "bo", markersize=1)
    !isinteractive() && show()
    sleep(0.2)
end

const d = BinaryCoding(6, [-3, -3], [3, 3]) #Two variables -3 <= x_i <= 3, encoded with a precision of 6 decimal places

z1(x1, x2) = -(3(1-x1)^2 * exp(-x1^2 - (x2+1)^2) - 10(x1/5 - x1^3 - x2^5) * exp(-x1^2-x2^2) -3exp(-(x1+2)^2 - x2^2) + 0.5(2x1 + x2))
z2(x1, x2) = -(3(1+x2)^2 * exp(-x2^2 - (1-x1)^2) - 10(-x2/5 + x2^3 + x1^5) * exp(-x1^2-x2^2) - 3exp(-(2-x2)^2 - x1^2))
z(x) = z1(x[1], x[2]), z2(x[1], x[2])

nsga(300, 50, z, d, fplot = plot_pop, seed = [[.0, 0.], [1., 1.5]])
```

* You don't have to provide a initialisation function anymore, a bitstring of the appropriate length will be generated.
* The seed can be passed as a vector of phenotypes, not a vector of genotypes

You can also use `BinaryCoding(ϵ::Int, types, lb, ub)` to encode a mix of integer, continuous or binary variables, with `types` a vector of symbols : `( :Int |  :Cont | :Bin )`, 
