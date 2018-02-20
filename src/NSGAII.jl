__precompile__()
module NSGAII

export nsga, nsga_max, MixedCoding, RealCoding
using ProgressMeter, StaticArrays, Compat, Compat.Random
using FunctionWrappers: FunctionWrapper


struct NSGA_FUNCS{G, P, N, Y}
    feval::FunctionWrapper{NTuple{N,Y}, Tuple{P}}
    finit::FunctionWrapper{G, Tuple{}}
    fCV::FunctionWrapper{Float64, Tuple{P}}
    fcross::FunctionWrapper{Nothing, Tuple{G,G,G,G}}
    fmut::FunctionWrapper{Nothing, Tuple{G}}
    fdecode::FunctionWrapper{P, Tuple{G}}

    function NSGA_FUNCS(::Type{G}, ::Type{P}, ::Val{N}, ::Type{Y}, feval, finit, fCV, fcross, fmut, fdecode) where {G,P,N,Y}
        new{G,P,N,Y}( FunctionWrapper{NTuple{N,Y}, Tuple{P}}(feval),
                    FunctionWrapper{G, Tuple{}}(finit),
                    FunctionWrapper{Float64, Tuple{P}}(fCV),
                    FunctionWrapper{Nothing, Tuple{G,G,G,G}}(fcross),
                    FunctionWrapper{Nothing, Tuple{G}}(fmut),
                    FunctionWrapper{P, Tuple{G}}(fdecode),
                    )
    end
end


include("indivs.jl")
include("functions.jl")
include("crossover.jl")
include("mutation.jl")
include("mixedcoding.jl")
include("vOptWrapper.jl")




function nsga(popSize::Integer, nbGen::Integer, z::Function, init::Function ; 
    fdecode=identity, fCV = x->0., pmut= 0.05, fmut=default_mutation!, 
    fcross = default_crossover!, seed=typeof(init())[], fplot = (x)->nothing, plotevery=10)
    g = init() ; G = typeof(g)
    p = fdecode(g) ; P = typeof(p)
    y = z(p) ; N = length(y) ; Y = eltype(y)
    funcs = NSGA_FUNCS(G,P,Val(N),Y, z,init,fCV,fcross,fmut,fdecode)
    return _nsga(funcs, Min(), popSize, nbGen, pmut, seed, fplot, plotevery)
end

function nsga(popSize::Integer, nbGen::Integer, z::Function, mc::MixedCoding, init::Function=()->bitrand(mc.nbbitstotal); 
    fCV = x->0., pmut= 0.05, fmut=default_mutation!, fcross = default_crossover!, 
    seed=Vector{Float64}[], fplot = (x)->nothing, plotevery=10)
    g = init() ; G = typeof(g)
    p = decode(g, mc) ; P = typeof(p)
    y = z(p) ; N = length(y) ; Y = eltype(y)
    funcs = NSGA_FUNCS(G,P,Val(N),Y, z,init,fCV,fcross,fmut,x->decode(x, mc))
    return _nsga(funcs, Min(), popSize, nbGen, pmut, encode.(seed, mc), fplot, plotevery)
end

function nsga_max(popSize::Integer, nbGen::Integer, z::Function, init::Function ; 
    fdecode=identity, fCV = x->0., pmut= 0.05, fmut=default_mutation!, 
    fcross = default_crossover!, seed=typeof(init())[], fplot = (x)->nothing, plotevery=10)
    g = init() ; G = typeof(g)
    p = fdecode(g) ; P = typeof(p)
    y = z(p) ; N = length(y) ; Y = eltype(y)
    funcs = NSGA_FUNCS(G,P,Val(N),Y, z,init,fCV,fcross,fmut,fdecode)
    return _nsga(funcs, Max(), popSize, nbGen, pmut, seed, fplot, plotevery)
end

function nsga_max(popSize::Integer, nbGen::Integer, z::Function, mc::MixedCoding, init::Function=()->bitrand(mc.nbbitstotal); 
    fCV = x->0., pmut= 0.05, fmut=default_mutation!, fcross = default_crossover!, 
    seed=Vector{Float64}[], fplot = (x)->nothing, plotevery=10)
    g = init() ; G = typeof(g)
    p = decode(g, mc) ; P = typeof(p)
    y = z(p) ; N = length(y) ; Y = eltype(y)
    funcs = NSGA_FUNCS(G,P,Val(N),Y, z,init,fCV,fcross,fmut,x->decode(x, mc))
    return _nsga(funcs, Max(), popSize, nbGen, pmut, encode.(seed, mc), fplot, plotevery)
end


end # module
