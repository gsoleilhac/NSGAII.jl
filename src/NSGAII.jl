__precompile__()
module NSGAII

export nsga, nsga_max, nsga_binary, MixedCoding, RealCoding
using ProgressMeter, StaticArrays, Compat, Compat.Random, Traceur

include("indivs.jl")
include("functions.jl")
include("crossover.jl")
include("mutation.jl")
include("mixedcoding.jl")
include("vOptWrapper.jl")

function nsga(popSize::Integer, nbGen::Integer, z::Function, init::Function ; 
    fdecode=identity, fdecode! = (geno,pheno)-> pheno.=geno, fCV = x->0., pmut= 0.05, fmut=default_mutation!, 
    fcross = default_crossover!, seed=typeof(init())[], fplot = (x)->nothing, plotevery=10)
	X = create_indiv(init(), fdecode, z, fCV)
    return _nsga(X, Min(), popSize, nbGen, init, z, fdecode, fdecode!, fCV , pmut, fmut, fcross, seed, fplot, plotevery)
end

function nsga(popSize::Integer, nbGen::Integer, z::Function, mc::MixedCoding, init::Function=()->bitrand(mc.nbbitstotal); 
    fCV = x->0., pmut= 0.05, fmut=default_mutation!, fcross = default_crossover!, 
    seed=Vector{Float64}[], fplot = (x)->nothing, plotevery=10)
    X = create_indiv(init(), x->decode(x, mc), z, fCV)
    return _nsga(X, Min(), popSize, nbGen, init, z, x->decode(x, mc), (g,f)->decode!(g, mc, f), fCV , pmut, fmut, fcross, encode.(seed, mc), fplot, plotevery)
end

function nsga_max(popSize::Integer, nbGen::Integer, z::Function, init::Function ; 
    fdecode=identity, fdecode! = (geno,pheno)-> pheno.=geno, fCV = x->0., pmut= 0.05, fmut=default_mutation!, 
    fcross = default_crossover!, seed=typeof(init())[], fplot = (x)->nothing, plotevery=10)
	X = create_indiv(init(), fdecode, z, fCV)
    return _nsga(X, Max(), popSize, nbGen, init, z, fdecode, fdecode!, fCV , pmut, fmut, fcross, seed, fplot, plotevery)
end

function nsga_max(popSize::Integer, nbGen::Integer, z::Function, mc::MixedCoding, init::Function=()->bitrand(mc.nbbitstotal); 
    fCV = x->0., pmut= 0.05, fmut=default_mutation!, fcross = default_crossover!, 
    seed=Vector{Float64}[], fplot = (x)->nothing, plotevery=10)
	X = create_indiv(init(), x->decode(x, mc), z, fCV)
    return _nsga(X, Max(), popSize, nbGen, init, z, x->decode(x, mc), (g,f)->decode!(g, mc, f), fCV , pmut, fmut, fcross, encode.(seed, mc), fplot, plotevery)
end

end # module
