module NSGAII
export nsga, MixedCoding, RealCoding
using ProgressMeter, StaticArrays, Compat, Compat.Random

if VERSION >= v"0.7-"
    Base.digits(T,n,base,pad) = digits(T, n, base=base, pad=pad)
end

include("indivs.jl")
include("functions.jl")
include("crossover.jl")
include("mutation.jl")
include("mixedcoding.jl")
include("vOptWrapper.jl")


function nsga(popSize::Integer, nbGen::Integer, z::Function, init::Function ; 
    fdecode=identity, fCV = x->0., pmut= 0.05, fmut=default_mutation!, 
    fcross = default_crossover, seed=typeof(init())[], fplot = (x)->nothing)
    return _nsga(popSize, nbGen, init, z, fdecode, fCV , pmut, fmut, fcross, seed, fplot)
end

function nsga(popSize::Integer, nbGen::Integer, z::Function, mc::MixedCoding, init::Function=()->bitrand(mc.nbbitstotal); 
    fCV = x->0., pmut= 0.05, fmut=default_mutation!, fcross = default_crossover, 
    seed=Vector{Float64}[], fplot = (x)->nothing)
    return _nsga(popSize, nbGen, init, z, x->decode(x, mc), fCV , pmut, fmut, fcross, encode.(seed, mc), fplot)
end

end # module
