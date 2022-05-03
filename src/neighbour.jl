#=
Code for building neighbour lists and extended point array
=#
import Base
export ExtendedPointArray, NeighbourList, eachneighbour, nions_extended, nions_orig, num_neighbours
export rebuild!, update!


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
    rcut::Float64
    lattice::Matrix{Float64}
end

function Base.show(io::IO, s::ExtendedPointArray)
    print(io, "ExtendedPointArray of $(length(s.indices)) points from $(length(s.orig_positions)) points")
end

"""
    ExtendedPointArray(cell::Cell, rcut)

Constructed an ExtendedPointArray from a given structure
"""
function ExtendedPointArray(cell::Cell, rcut)
    rcut = convert(Float64, rcut)
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
    ExtendedPointArray(indices, shiftidx, shifts, [SVector{3}(x) for x in eachcol(pos_extended)], original_positions, 
                       inoshift, rcut, copy(cellmat(lattice(cell))))
end

"""
Rebuild the ExtendedPointArray for an existing cell
"""
function rebuild!(ea::ExtendedPointArray, cell)
    lattice_change = !all(ea.lattice .== cellmat(cell))
    i = 1 
    # Rebuild shift vectors from scratch is the lattice vectors have been changed
    if lattice_change
        newshifts = CellBase.shift_vectors(cellmat(lattice(cell)), ea.rcut;safe=false)
        # Number of vectors change
        nnew = length(newshifts)
        nold = length(ea.shiftvecs)
        dl = nnew - nold
        if dl != 0
            resize!(ea.shiftvecs, nnew)
            # Also resize the positions array and indices
            ntot = nnew * length(ea.orig_positions)
            resize!(ea.positions, ntot)
            resize!(ea.indices, ntot)
            resize!(ea.shiftidx, ntot)
        end
        ea.shiftvecs .= newshifts

        # Rebuild the extended arrays
        ea.orig_positions .= sposarray(cell)
        for (idx, pos_orig) in enumerate(ea.orig_positions)   # Each original positions
            for (ishift, shiftvec) in enumerate(ea.shiftvecs)   # Each shift positions
                ea.positions[i] = pos_orig .+ shiftvec
                ea.indices[i] = idx
                ea.shiftidx[i] = ishift
                i += 1
            end
        end
    else
        # no change in lattice shift all extended points with the displacements of the original ones
        spos = sposarray(cell)
        m = 1
        for (ipos, pos) in enumerate(spos)
            disp = pos - ea.orig_positions[ipos]
            # Displace all images by this amount
            for _ = 1:length(ea.shiftvecs)
                ea.positions[m] = ea.positions[m] + disp
                m += 1
            end
        end
        ea.orig_positions .= spos
    end
    ea
end


"""
Type for representing a neighbour list
"""
struct NeighbourList{T, N}
    ea::ExtendedPointArray{T}
    "Extended indice of the neighbours"
    extended_indices::Matrix{Int}
    "Original indice of the neighbours"
    orig_indices::Matrix{Int}
    "Distance to the neighbours"
    distance::Matrix{Float64}
    "Vector displacement to the neighbours"
    vectors::Array{SVector{N, Float64}, 2}
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
    println(io, "NeighbourList of $(nions_orig(n)) ($(nions_extended(n)) extended)")
    print(io, "Current max neighbours $(maximum(n.nneigh))\n")
    print(io, "Cut off radius: $(n.rcut)")
end

"""
    NeighbourList(ea::ExtendedPointArray, rcut, nmax=100; savevec=false)

Construct a NeighbourList from an extended point array for the points in the original cell
"""
function NeighbourList(ea::ExtendedPointArray, rcut, nmax=100; savevec=false, ndim=3)
    rcut = convert(Float64, rcut) 
    norig = length(ea.orig_positions)
    extended_indices = zeros(Int, nmax, norig)
    orig_indices = zeros(Int, nmax, norig)
    distance = fill(-1., nmax, norig)
    nneigh = zeros(Int, norig)

    @assert rcut >= ea.rcut "Cut off radius is large that that of the periodic image cut off."

    # Save vectors or not
    base = @SVector fill(-1., ndim)
    savevec ? vectors = fill(base, nmax, norig) : vectors = fill(SA[-1., -1. , -1.], 1, 1)
    nl = NeighbourList(
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
    rebuild!(nl, ea)
    nl
end

"""
    rebuild!(nl::NeighbourList, ea::ExtendedPointArray)

Perform a full rebuild of the neighbour list from scratch
"""
function rebuild!(nl::NeighbourList, ea::ExtendedPointArray)
    
    extended_indices = nl.extended_indices
    orig_indices = nl.orig_indices
    distance = nl.distance
    nneigh = nl.nneigh 
    nmax = nl.nmax
    savevec = nl.has_vectors
    vectors = nl.vectors
    rcut = nl.rcut
    # Reset
    fill!(vectors, vectors[1] .* 0.)
    fill!(distance, -1.)
    fill!(orig_indices, 0)
    fill!(extended_indices, 0)
    fill!(nneigh, 0)

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
                    savevec && (vectors[ineigh, iorig] = posj .- posi)
                end
            end
        end
        # Store the total number of neighbours for this point
        if ineigh > nmax
            throw(ErrorException("Too many neighbours, please increase the value of nmax to at least $(ineigh)"))
        end
        nneigh[iorig] = ineigh
    end
    nl
end

"""
    rebuild(nl::NeighbourList, cell::Cell) 

Perform a full rebuild of the NeighbourList with the latest geometry of the cell
"""
function rebuild!(nl::NeighbourList, cell::Cell) 
    rebuild!(nl.ea, cell)
    rebuild!(nl, nl.ea)
end

"""
Update the calculated distances/vectors but do not rebuild the whole neighbour list
"""
function update!(nl::NeighbourList)
    ea = nl.ea
    for (isite, nneight) in enumerate(nl.nneigh)
        inn = 1
        posi = ea.orig_positions[isite]
        while inn <= nneight
            jid_ext = nl.extended_indices[inn, isite]
            posj = ea.positions[jid_ext]
            dij = distance_between(posj, posi)
            # Update distance
            nl.distance[inn, isite] = dij
            # Update the vector
            nl.has_vectors && (nl.vectors[inn, isite] = posj .- posi)
            inn += 1
        end
    end
end

"""
    update!(nl::NeighbourList, cell::Cell) 

Update the NeighbourList with the latest geometry of the Cell.
No rebuilding is performed
"""
function update!(nl::NeighbourList, cell::Cell) 
    rebuild!(nl.ea, cell)
    update!(nl)
end


NeighbourList(cell::Cell, rcut, nmax=100;savevec=false) = NeighbourList(ExtendedPointArray(cell, rcut), rcut, nmax;savevec)

"Number of neighbours for a point"
num_neighbours(nl::NeighbourList, iorig) = nl.nneigh[iorig]


abstract type AbstractNLIterator end

"Iterator interface for going through all neighbours"
struct NLIterator{T} <: AbstractNLIterator
    nl::T
    iorig::Int
end

"Iterator interface for going through all neighbours"
struct NLIteratorWithVector{T} <: AbstractNLIterator
    nl::T
    iorig::Int
end


Base.length(nli::AbstractNLIterator) =  num_neighbours(nli.nl, nli.iorig)

function Base.iterate(nli::NLIterator, state=1)
    nl = nli.nl
    iorig = nli.iorig
    if state > nl.nneigh[nli.iorig] 
        return nothing
    end
    return (nl.orig_indices[state, iorig], nl.extended_indices[state, iorig], nl.distance[state, iorig]), state + 1
end

function Base.iterate(nli::NLIteratorWithVector, state=1)
    nl = nli.nl
    iorig = nli.iorig
    if state > nl.nneigh[nli.iorig] 
        return nothing
    end
    return (nl.orig_indices[state, iorig], nl.extended_indices[state, iorig], nl.distance[state, iorig], nl.vectors[state, iorig]), state + 1
end

"""
Iterate the neighbours of a site in the original cell.
Returns a tuple of (original_index, extended_index, distance) for each iteration
"""
eachneighbour(nl::NeighbourList, iorig) = NLIterator(nl, iorig)


"""
Iterate the neighbours of a site in the original cell.
Returns a tuple of (original_index, extended_index, distance, vector) for each iteration
"""
function eachneighbourvector(nl::NeighbourList, iorig) 
    @assert nl.has_vectors "NeighbourList is not build with distance vectors"
    NLIteratorWithVector(nl, iorig)
end