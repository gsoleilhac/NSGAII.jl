__precompile__()
module NSGAII

export nsga, nsga_max, nsga_binary, BinaryCoding
using ProgressMeter, Compat, Compat.Random

include("indivs.jl")
include("functions.jl")
include("crossover.jl")
include("mutation.jl")
include("binarycoding.jl")
include("vOptWrapper.jl")

function nsga(popSize::Integer, nbGen::Integer, z::Function, init::Function ; 
    fdecode=identity, fdecode! = (geno,pheno)-> pheno.=geno, fCV = x->0., pmut= 0.05, fmut=default_mutation!, 
    fcross = default_crossover!, seed=typeof(init())[], fplot = (x)->nothing, plotevery=1)
	X = create_indiv(init(), fdecode, z, fCV)
    return _nsga(X, Min(), popSize, nbGen, init, z, fdecode, fdecode!, fCV , pmut, fmut, fcross, seed, fplot, plotevery)
end

function nsga(popSize::Integer, nbGen::Integer, z::Function, bc::BinaryCoding, init::Function=()->bitrand(bc.nbbitstotal); 
    fCV = x->0., pmut= 0.05, fmut=default_mutation!, fcross = default_crossover!, 
    seed=Vector{Float64}[], fplot = (x)->nothing, plotevery=1)
    X = create_indiv(init(), x->decode(x, bc), z, fCV)
    return _nsga(X, Min(), popSize, nbGen, init, z, x->decode(x, bc), (g,f)->decode!(g, bc, f), fCV , pmut, fmut, fcross, encode.(seed, bc), fplot, plotevery)
end

function nsga_max(popSize::Integer, nbGen::Integer, z::Function, init::Function ; 
    fdecode=identity, fdecode! = (geno,pheno)-> pheno.=geno, fCV = x->0., pmut= 0.05, fmut=default_mutation!, 
    fcross = default_crossover!, seed=typeof(init())[], fplot = (x)->nothing, plotevery=1)
	X = create_indiv(init(), fdecode, z, fCV)
    return _nsga(X, Max(), popSize, nbGen, init, z, fdecode, fdecode!, fCV , pmut, fmut, fcross, seed, fplot, plotevery)
end

function nsga_max(popSize::Integer, nbGen::Integer, z::Function, bc::BinaryCoding, init::Function=()->bitrand(bc.nbbitstotal); 
    fCV = x->0., pmut= 0.05, fmut=default_mutation!, fcross = default_crossover!, 
    seed=Vector{Float64}[], fplot = (x)->nothing, plotevery=1)
	X = create_indiv(init(), x->decode(x, bc), z, fCV)
    return _nsga(X, Max(), popSize, nbGen, init, z, x->decode(x, bc), (g,f)->decode!(g, bc, f), fCV , pmut, fmut, fcross, encode.(seed, bc), fplot, plotevery)
end

end # module
