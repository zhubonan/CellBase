#=
For handling compositions
=#
import Base

"""
Type representing a composition
"""
struct Composition
    species::Vector{Symbol}
    counts::Vector{Float64}
    function Composition(species, counts)
        length(species) != length(counts) && error("Size mismatch") 
        sidx = sortperm(species)
        new(species[sidx], counts[sidx])
    end
end

"""
    Composition(pairs::Pair{Symbol, T}...) where T

Return `Composition` from pairs.
"""
function Composition(pairs::Pair{Symbol, T}...) where T
    x = Array{Symbol}(undef, length(pairs))
    y = Array{Float64}(undef, length(pairs))
    i = 1
    for (a, b) in pairs
        x[i] = a
        y[i] = b
        i += 1
    end
    Composition(x, y)
end

"""
    Composition(comp::Dict)

Return `Composition` from dictionary
"""
function Composition(comp::Dict)
    Composition(collect(comp)...)
end


function Base.getindex(t::Composition, key) 
    idx = findfirst(x -> x == key, t.species)
    isnothing(idx) && return 0.
    t.counts[idx]
end

Base.haskey(t::Composition, key) = key in t.species

Base.keys(t::Composition) = t.species
Base.pairs(t::Composition) = [a => b for (a, b) in zip(t.species, t.counts)]
Base.iterate(t::Composition) = Base.iterate(Base.pairs(t))
Base.iterate(t::Composition, i) = Base.iterate(Base.pairs(t), i)
Base.length(t::Composition) = Base.length(t.species)

"""
    Composition(cell::Cell)

Return the composition of a Cell
"""
function Composition(cell::Cell)
    sp_array = species(cell)
    unique_sp = sort(unique(sp_array))
    num_atoms = Array{Int}(undef, size(unique_sp))
    for i in 1:length(unique_sp)
       num_atoms[i] = count(x -> x == unique_sp[i], sp_array)
    end
    Composition(unique_sp, num_atoms)
end