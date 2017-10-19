# NSGAII

[![Build Status](https://travis-ci.org/gsoleilhac/NSGAII.jl.svg?branch=master)](https://travis-ci.org/gsoleilhac/NSGAII.jl)

[![codecov.io](http://codecov.io/github/gsoleilhac/NSGAII.jl/coverage.svg?branch=master)](http://codecov.io/github/gsoleilhac/NSGAII.jl?branch=master)

*A generic NSGAII implementation in Julia*

## Installation
`Pkg.clone("https://github.com/gsoleilhac/NSGAII.jl")`


## Usage : 

```
popsize = 200
nbGenerations = 100
pMutation = 0.05

#define how to generate a random genotype : 
init_function = () -> randperm(N) #e.g : a permutation coding

#define how to evaluate a genotype : 
z(x) = z1(x), z2(x) ... # Must return a Tuple

nsga(popsize, nbGenerations, init_function, eval, pMutation)
```

By default, a PMX Crossover and a random swap mutation will be used if the genotype is a Vector{Int}

and a two-point crossover and random flips will be used for Vector{Bool} / BitArrays.

You can define your own crossover / mutations operators with the keywords `fmut` and `fcross`

```
nsga(popsize, nbGen, init, z , fmut = ..., fcross = ...)
```

You can also define some starting solutions with the keyword `seed`
e.g : 
```
x1, x2, x3 = greedy(...)
nsga(... , seed = [x1, x2, x3])
```

Finally a function to plot the population during the evolution can be passed with the keyword `fplot` ; 
the genotype of an individuam is accessed with `indiv.x` and it's objective function values with `indiv.y`


## RealCoding

The package also exports `RealCoding(eps, lb, ub)` and `decode(x, d::RealCoding)` which can be used to easily represent Real values with eps decimal places.

`encode(x, d)` is also provided for the seed population

```
d = RealCoding(6, [-3, -3], [6, 6]) # Codes two reals with 6 decimal places with lower bound -3 and upper bound 6
z1(x1, x2) = -(3(1-x1)^2 * exp(-x1^2 - (x2+1)^2) - ...)
z2(x1, x2) = -(3(1+x2)^2 * exp(-x2^2 - (1-x1)^2) - ...)
z(x) = begin 
    x1, x2 = decode(x, d)
    z1(x1, x2), z2(x1, x2)
end
nsga(100, 50, ()->bitrand(d.nbbitstotal), z)
```
