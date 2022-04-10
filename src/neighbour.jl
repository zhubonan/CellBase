#=
Code for building neighbour lists and extended point array
=#
import Base
export ExtendedPointArray, NeighbourList, eachneighbour, nions_extended, nions_orig, num_neighbours


"""
Represent an array of points after expansion by periodic boundary
"""
struct ExtendedPointArray{T}
    "Original point indices"
    indices::Vector{Int}
    "Index of the shift"
    shiftidx::Vector{Int}
    "Shift vectors"
    shiftvecs::Vector{SVector{3, Float64}}
    "Positions"
    positions::Vector{T}
    "original Positions"
    orig_positions::Vector{T}
    "Index of of the all-zero shift vector"
    inoshift::Int
end

function Base.show(io::IO, s::ExtendedPointArray)
    print(io, "ExtendedPointArray of $(length(s.indices)) points from $(length(s.orig_positions)) points")
end

"""
    ExtendedPointArray(cell::Cell, rcut)

Constructed an ExtendedPointArray from a given structure
"""
function ExtendedPointArray(cell::Cell, rcut)

    ni = nions(cell)
    shifts = CellBase.shift_vectors(cellmat(lattice(cell)), rcut;safe=false)
    indices = zeros(Int, ni * length(shifts))
    shiftidx = zeros(Int, ni * length(shifts))
    pos_extended = zeros(eltype(positions(cell)), 3, ni * length(shifts))
    inoshift = findfirst(x -> all( x .== 0.), shifts)

    i = 1
    original_positions = sposarray(cell)
    for (idx, pos_orig) in enumerate(original_positions)   # Each original positions
        for (ishift, shiftvec) in enumerate(shifts)   # Each shift positions
            pos_extended[:, i] .= pos_orig .+ shiftvec
            indices[i] = idx
            shiftidx[i] = ishift
            i += 1
        end
    end
    ExtendedPointArray(indices, shiftidx, shifts, [SVector{3}(x) for x in eachcol(pos_extended)], original_positions, inoshift)
end


"""
Type for representing a neighbour list
"""
struct NeighbourList{T}
    ea::ExtendedPointArray{T}
    "Extended indice of the neighbours"
    extended_indices::Matrix{Int}
    "Original indice of the neighbours"
    orig_indices::Matrix{Int}
    "Distance to the neighbours"
    distance::Matrix{Float64}
    "Vector displacement to the neighbours"
    vectors::Array{Float64, 3}
    "Number of neighbours"
    nneigh::Vector{Int}
    "Maximum number of neighbours that can be stored"
    nmax::Int
    "Contains vector displacements or not"
    has_vectors::Bool
    rcut::Float64
end

"Number of ions in the original cell"
nions_orig(n::ExtendedPointArray) = length(n.orig_positions)
"Number of ions in the original cell"
nions_orig(n::NeighbourList) = nions_orig(n.ea) 

"Number of ions in the extended cell"
nions_extended(n::ExtendedPointArray) = length(n.positions)
"Number of ions in the extended cell"
nions_extended(n::NeighbourList) = nions_extended(n.ea)


function Base.show(io::IO, n::NeighbourList)
    print("NeighbourList of maximum sizs $(n.nmax) for $(nions_orig(n))/$(nions_extended(n)) atoms")
end

"""
    NeighbourList(ea::ExtendedPointArray, rcut, nmax=100; savevec=false)

Construct a NeighbourList from an extended point array for the points in the original cell
"""
function NeighbourList(ea::ExtendedPointArray, rcut, nmax=100; savevec=false)
    
    norig = length(ea.orig_positions)
    extended_indices = zeros(Int, nmax, norig)
    orig_indices = zeros(Int, nmax, norig)
    distance = fill(-1., nmax, norig)
    nneigh = zeros(Int, norig)

    # Save vectors or not
    savevec ? vectors = fill(-1., 3, norig, nmax) : vectors = fill(-1., 1, 1, 1)

    for (iorig, posi) in enumerate(ea.orig_positions)
        ineigh = 0
        for (j, posj) in enumerate(ea.positions)
            (ea.indices[j] == iorig) && (ea.shiftidx[j] == ea.inoshift) && continue
            dist = distance_between(posj, posi)
            # Store the information if the distance is smaller than the cut off
            if dist < rcut
                ineigh += 1
                if ineigh <= nmax
                    distance[ineigh, iorig] = dist
                    # Store the extended index
                    extended_indices[ineigh, iorig] = j
                    # Store the index of the point in the original cell
                    orig_indices[ineigh, iorig] = ea.indices[j]
                    savevec && (vectors[:, ineigh, iorig] .= posj .- posi)
                end
            end
        end
        # Store the total number of neighbours for this point
        if ineigh > nmax
            throw(ErrorException("Too many neighbours, please increase the value of nmax to at least $(ineigh)"))
        end
        nneigh[iorig] = ineigh
    end
    NeighbourList(
        ea,
        extended_indices,
        orig_indices,
        distance,
        vectors,
        nneigh,
        nmax,
        savevec,
        rcut,
    )
end


NeighbourList(cell::Cell, rcut, nmax=100;savevec=false) = NeighbourList(ExtendedPointArray(cell, rcut), rcut, nmax;savevec)

"Number of neighbours for a point"
num_neighbours(nl::NeighbourList, iorig) = nl.nneigh[iorig]


"Iterator interface for going through all neighbours"
struct NLIterator
    nl::NeighbourList
    iorig::Int
end

Base.length(nli::NLIterator) =  num_neighbours(nli.nl, nli.iorig)

function Base.iterate(nli::NLIterator, state=1)
    nl = nli.nl
    iorig = nli.iorig
    if state > nl.nneigh[nli.iorig] 
        return nothing
    end
    return (nl.orig_indices[state, iorig], nl.extended_indices[state, iorig], nl.distance[state, iorig]), state + 1
end

"""
Iterate the neighbours of a site in the original cell.
Returns a tuple of (original_index, extended_index, distance) for each iteration
"""
eachneighbour(nl::NeighbourList, iorig) = NLIterator(nl, iorig)