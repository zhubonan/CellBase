#=
Representation of single atomic Sites

The `Site` type is deprecated.
=#

using Printf
using StaticArrays
export Site, distance_between, distance_squared_between

"""
Site{T}

Represent a site with absolute coordinates
"""
struct Site{T}
    position::T
    index::Int
    symbol::Symbol
end

"Initialise a site from only position and symbol"
function Site(pos, symbol)
    Site(pos, 0, symbol)
end

function Base.show(io::IO, s::Site{T}) where T
    print(io, @sprintf "Site{%s}: %5s %4d at %10.5f %10.5f %10.5f" T s.symbol s.index s.position[1] s.position[2] s.position[3])
end

function Base.show(io::IO, ::MIME"text/plain", s::Site{T}) where T
    print(io, @sprintf "Site{%s}:\n" T)
    print(io, @sprintf "%5s %4d %10.5f %10.5f %10.5f" s.symbol s.index s.position[1] s.position[2] s.position[3])
end

Base.position(s::Site) = s.position
index(s::Site) = s.index
symbol(s::Site) = s.symbol


distance_squared_between(s1::Site, s2::Site) = sum((s1.position .- s2.position) .^ 2)
distance_squared_between(s1::AbstractVector, s2::AbstractVector) = sum((s1 - s2) .^ 2) 
distance_squared_between(s1::AbstractVector, s2::AbstractVector, shift) = sum((s1 - s2 .- shift) .^ 2) 
distance_squared_between(s1::Site, s2::Site, shift) = sum((s1.position .- s2.position .- shift) .^ 2)

distance_between(s1, s2, shift) = sqrt(distance_squared_between(s1, s2, shift))
distance_between(s1, s2) = sqrt(distance_squared_between(s1, s2))


"Unit vector from site 1 to shifted site 2"
function unit_vector_between(s1::Site, s2::Site, shift)
    vtmp .= s2.position .- s1.position .+ shift
    vtmp ./ norm(vtmp)
end