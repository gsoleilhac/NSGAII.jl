__precompile__()
module NSGAII

export nsga, nsga_max, MixedCoding, RealCoding
using ProgressMeter, StaticArrays, Compat, Compat.Random

include("indivs.jl")
include("functions.jl")
include("crossover.jl")
include("mutation.jl")
include("mixedcoding.jl")
include("vOptWrapper.jl")

function nsga(popSize::Integer, nbGen::Integer, z::Function, init::Function ; 
    fdecode=identity, fCV = x->0., pmut= 0.05, fmut=default_mutation!, 
    fcross = default_crossover!, seed=typeof(init())[], fplot = (x)->nothing, plotevery=10)
	X = typeof(indiv(init(), fdecode, z, fCV))
    return _nsga(X, Min(), popSize, nbGen, init, z, fdecode, fCV , pmut, fmut, fcross, seed, fplot, plotevery)
end

function nsga(popSize::Integer, nbGen::Integer, z::Function, mc::MixedCoding, init::Function=()->bitrand(mc.nbbitstotal); 
    fCV = x->0., pmut= 0.05, fmut=default_mutation!, fcross = default_crossover!, 
    seed=Vector{Float64}[], fplot = (x)->nothing, plotevery=10)
	X = typeof(indiv(init(), x->decode(x, mc), z, fCV))
    return _nsga(X, Min(), popSize, nbGen, init, z, x->decode(x, mc), fCV , pmut, fmut, fcross, encode.(seed, mc), fplot, plotevery)
end

function nsga_max(popSize::Integer, nbGen::Integer, z::Function, init::Function ; 
    fdecode=identity, fCV = x->0., pmut= 0.05, fmut=default_mutation!, 
    fcross = default_crossover!, seed=typeof(init())[], fplot = (x)->nothing, plotevery=10)
	X = typeof(indiv(init(), fdecode, z, fCV))
    return _nsga(X, Max(), popSize, nbGen, init, z, fdecode, fCV , pmut, fmut, fcross, seed, fplot, plotevery)
end

function nsga_max(popSize::Integer, nbGen::Integer, z::Function, mc::MixedCoding, init::Function=()->bitrand(mc.nbbitstotal); 
    fCV = x->0., pmut= 0.05, fmut=default_mutation!, fcross = default_crossover!, 
    seed=Vector{Float64}[], fplot = (x)->nothing, plotevery=10)
	X = typeof(indiv(init(), x->decode(x, mc), z, fCV))
    return _nsga(X, Max(), popSize, nbGen, init, z, x->decode(x, mc), fCV , pmut, fmut, fcross, encode.(seed, mc), fplot, plotevery)
end

end # module
