__precompile__()
module NSGAII

export nsga, nsga_max, BinaryCoding

using ProgressMeter
using Random
using LinearAlgebra

include("indivs.jl")
include("functions.jl")
include("crossover.jl")
include("mutation.jl")
include("binarycoding.jl")
# include("vOptWrapper.jl")

function nsga(popSize::Integer, nbGen::Integer, z::Function, init::Function ; 
    fCV = x -> 0., pmut = 0.05, fmut = default_mutation!, fcross = default_crossover!, 
    seed = typeof(init())[], fplot = x -> nothing, plotevery = 1, showprogress = true)

    X = createIndiv(init(), z, fCV)
    
    return _nsga(X, Min(), popSize, nbGen, init, z, fCV, pmut, fmut, fcross, seed, fplot, plotevery, showprogress ? 0.5 : Inf)
end

function nsga(popSize::Integer, nbGen::Integer, z::Function, bc::BinaryCoding ; 
    fCV = x -> 0., pmut = 0.05, fmut = indiv -> default_mutation!(indiv.x), fcross = (indivs...) -> default_crossover!(getproperty.(indivs, :x)...), 
    seed = Vector{Float64}[], fplot = x -> nothing, plotevery = 1, showprogress = true)

    init = () -> BinaryCodedIndiv(bitrand(bc.nbbitstotal), zeros(bc.nbvar))
    _z = indiv -> (decode!(indiv, bc) ; z(indiv.p))

    X = createIndiv(init(), _z, indiv -> fCV(indiv.p))
    return _nsga(X, Min(), popSize, nbGen, init, _z, fCV , pmut, fmut, fcross, encode.(seed, Ref(bc)), fplot, plotevery, showprogress ? 0.5 : Inf)
end

function nsga_max(popSize::Integer, nbGen::Integer, z::Function, init::Function ; 
    fCV = x -> 0., pmut = 0.05, fmut = default_mutation!, fcross = default_crossover!, 
    seed = typeof(init())[], fplot = x -> nothing, plotevery = 1, showprogress = true)

    X = createIndiv(init(), z, fCV)
    
    return _nsga(X, Max(), popSize, nbGen, init, z, fCV, pmut, fmut, fcross, seed, fplot, plotevery, showprogress ? 0.5 : Inf)
end

function nsga_max(popSize::Integer, nbGen::Integer, z::Function, bc::BinaryCoding ; 
    fCV = x -> 0., pmut = 0.05, fmut = indiv -> default_mutation!(indiv.x), fcross = (indivs...) -> default_crossover!(getproperty.(indivs, :x)...), 
    seed = Vector{Float64}[], fplot = x -> nothing, plotevery = 1, showprogress = true)

    init = () -> BinaryCodedIndiv(bitrand(bc.nbbitstotal), zeros(bc.nbvar))
    _z = indiv -> (decode!(indiv, bc) ; z(indiv.p))

    X = createIndiv(init(), _z, indiv -> fCV(indiv.p))
    return _nsga(X, Max(), popSize, nbGen, init, _z, fCV , pmut, fmut, fcross, encode.(seed, Ref(bc)), fplot, plotevery, showprogress ? 0.5 : Inf)
end


end # module
