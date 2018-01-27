# NSGAII

[![Build Status](https://travis-ci.org/gsoleilhac/NSGAII.jl.svg?branch=master)](https://travis-ci.org/gsoleilhac/NSGAII.jl)

[![codecov.io](http://codecov.io/github/gsoleilhac/NSGAII.jl/coverage.svg?branch=master)](http://codecov.io/github/gsoleilhac/NSGAII.jl?branch=master)

## Installation
`Pkg.clone("https://github.com/gsoleilhac/NSGAII.jl")`


## Usage : 

```
popsize = 200
nbGenerations = 100

#define how to generate a random genotype : 
init() = randperm(N) #e.g : a permutation coding

#define how to evaluate a phenotype :
z(x) = z1(x), z2(x) ... # Must return a Tuple

nsga(popsize, nbGenerations, z, init)

#By default, the identity function is used to calculate the phenotype from the genotype but this can be changed with the keyword fdecode

#A constraint violation function can be passed with the keyword fCV
#It should return 0. if the solution is feasible and a value > 0 otherwise.

#Mutation probability can be changed with the keyword pmut (default is 5%)

nsga(popsize, nbGen, z, init, fCV = (genotype) -> ... , fdecode = (phenotype) -> ... pmut = 0.2)

```

By default, a PMX Crossover and a random swap mutation will be used if the genotype is a Vector{Int}

and a two-point crossover and random flips will be used for Vector{Bool} / BitArrays.

Other crossover / mutations operators can be passed with the keywords `fmut` and `fcross`

```
nsga(popsize, nbGen, z, init, fmut = ..., fcross = ...)
```

Starting solutions are defined with the keyword `seed` ; they must be valid genotypes
e.g : 
```
x1, x2, x3 = greedy(...)
nsga(popsize, nbGen, z, init, seed = [x1, x2, x3])
```

A plot function can be passed with the keyword `fplot`
e.g. with PyPlot:
```
function plot_pop(P)
    clf()
    p = plot(map(x -> x.y[1], P), map(x -> x.y[2], P), "bo", markersize=1)
    sleep(0.2)
end

nsga(popsize, nbGen, z, init, fplot = plot_pop)
```

## RealCoding

`RealCoding(eps, lb, ub)` can be used to easily represent Real values with eps decimal places.

```
rc = RealCoding(5, [-4, -4], [6, 6]) #Codes two reals -4 <= x <= 6 with a precision of 5 decimal places
nsga(100, 50, z, rc, seed=[(-3., 6), (2.5, 0)])
#Note that the initilisation function is not needed anymore, a random sequence of bits with the appropriate length will be generated
#although it can still be overwritten : nsga(100, 50, z, rc, init, ...)

#The seed passed must be a phenotype here, it will automatically be encoded
```


## vOptGeneric

This package supports models declared with JuMP / vOptGeneric, with the restriction that all variables must be bounded.

See examples/Mavrota_MILP.jl
