using Printf

export Site, distance_between, distance_squared_between, displace!
"""
Site{T}

Represent a site with absolute coordinates
"""
struct Site{T}
    position::Vector{T}
    index::Int
    symbol::Symbol
end

"Initialise a site from only position and symbol"
function Site(pos::Vector{T}, symbol::Symbol) where T
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

"Distance between two sites"
function distance_between(s1::Site{T}, s2::Site{T}) where T
    sqrt(distance_squared_between(s1, s2))
end

"Distance between two sites with shifts applied to s2"
function distance_between(s1::Site{T}, s2::Site{T}, shift::Vector) where T
    sqrt(distance_squared_between(s1, s2, shift))
end

"Distance between two sites"
function distance_squared_between(s1::Site{T}, s2::Site{T}, shift::Vector) where T
    s = 0.0
    for i in 1:3
        s += (s1.position[i] - s2.position[i] - shift[i]) ^ 2
    end
    s
end

"Distance between two sites squared"
function distance_squared_between(s1::Site{T}, s2::Site{T}) where T
    s = 0.0
    for i in 1:3
        s += (s1.position[i] - s2.position[i]) ^ 2
    end
    s
end

"Distance between two three vectors"
function distance_squared_between(s1::AbstractVector, s2::AbstractVector) where T
    s = 0.0
    nv = length(s1)
    for i in 1:nv
        s += (s1[i] - s2[i]) ^ 2
    end
    s
end


"Unit vector from site 1 to shifted site 2"
function unit_vector_between!(vtmp::Vector{T}, s1::Site{T}, s2::Site{T}, shift::Vector) where T
    for i in 1:3
        vtmp[i] = s2.position[i] + shift[i] - s1.position[i] 
    end
    vtmp[:] /= norm(vtmp)
end

displace!(s::Site, vec::AbstractVector) = s.position[:] += vec
