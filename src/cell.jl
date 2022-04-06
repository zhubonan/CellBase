using Printf
using PeriodicTable
using LinearAlgebra

export Cell, nions, positions, species, lattice, volume, cellmat, cellvecs, wrap!, cellpar
"""
A Cell represents a periodic structure in three-dimensional space
"""
mutable struct Cell{T}
    lattice::Lattice{T}                 # Lattice of the structure
    sites::Vector{Site{T}}
    set_indices::Vector{Int}             # Masking array for the sets
    arrays::Dict{Symbol, AbstractArray}         # Any additional arrays
    info::Dict
end

"""
    Cell(l::Lattice, s::Vector{Site{T}}) where T

Construct a Cell type from lattice and a list of sites
"""
function Cell(l::Lattice, s::Vector{Site{T}}) where T
    set_indices = collect(1:length(s))
    arrays = Dict{Symbol, AbstractArray}()
    Cell(l, s, set_indices, arrays, Dict())
end

"""
    Cell(;lat::Lattice, pos::Matrix, numbers::Vector{Int})

Constructure the Cell type from lattice, positions and numbers
"""
function Cell(lat::Lattice, pos::Matrix{T}, numbers::Vector{Int}) where T
    @assert size(pos)[1] == 3 "Only support three dimensional structure"
    species = [Symbol(elements[i].symbol) for i in numbers]
    Cell(lat, pos, species)
end

"""
    Cell(;lat::Lattice, pos::Matrix, species::Vector{Symbol})

Constructure the Cell type from lattice, positions and species
"""
function Cell(lat::Lattice, pos::Matrix{T}, species::Vector{Symbol}) where T
    set_indices = [i for i in 1:size(pos)[2]]
    arrays = Dict{Symbol, AbstractArray}()
    sites = Site{T}[]
    for n in 1:size(pos)[2]
        push!(sites, Site(pos[:, n], n, species[n]))
    end
    Cell(lat, sites, set_indices, arrays, Dict())
end

"""
Clip a structure with a given indexing array
"""
function clip(s::Cell, mask::AbstractVector)
    new_sites = sites(s)[mask]
    new_sets = set_indices(s)[mask]  # Keep the original set index
    Cell(lattice(s), new_sites, new_sets, Dict{Symbol, AbstractArray}(), s.info)
end

# Basic interface 
"Number of atoms in a structure"
nions(structure::Cell) = length(structure.sites)

"Number of atoms in a structure"
natoms(structure::Cell) = length(structure.sites)

"Positions of ions in a structure"
function positions(structure::Cell{T}) where T
    posarray = Array{T}(undef, size(structure.sites[1].position)[1], nions(structure))
    for (i, site) in enumerate(structure.sites)
        posarray[:, i] = site.position
    end
    posarray
end

"Species names"
species(structure::Cell) = Symbol[site.symbol for site in structure.sites]

"Return the atomic numbers"
atomic_numbers(structure::Cell) = Int[elements[site.symbol].number for site in structure.sites]

## Wrapper for the Lattice ###

"Lattice type of a structure"
lattice(structure::Cell) = structure.lattice

"Cell parameters type of a structure"
cellpar(structure::Cell) = cellpar(lattice(structure))

"Volume of a structure"
volume(structure::Cell) = volume(lattice(structure))

"""Get the cell matrix of a structure with column vectors"""
cellmat(structure::Cell) = cellmat(lattice(structure))

## END ##


"Get additional array"
array(structure::Cell, arrayname::Symbol) = structure.arrays[arrayname]

"Names of the arrays in a Cell"
arraynames(structure::Cell) = keys(structure.arrays)

"Number of sets in a structures"
nsets(structure::Cell) = length(unique(structure.set_indices))

"""Indices for the sets"""
set_indices(structure::Cell) = structure.set_indices

"""Info dictionary"""
info(structure::Cell) = structure.info

"""Attach info dictionary"""
attachinfo!(structure::Cell, info::Dict) = structure.info = info

"""Number of formula units"""
num_fu(structure::Cell) = formula_and_factor(structure)[2]

"Reduced formula of a structures"
reduced_fu(structure::Cell) = formula_and_factor(structure)[1]

"""Sites of a structure"""
sites(structure::Cell) = structure.sites

"""
Sort symbols by atomic numbers
NOTE: dose not handle cases where the symbol is not an element!!!
"""
function sorted_symbols(symbols)
    z_array = [elements[sp].number for sp in symbols]
    perm = sortperm(z_array)
    return symbols[perm]
end

"""Reduced formula and number of units in a structure"""
function formula_and_factor(structure::Cell)

    # Check if computed results already exists
    info_dict = info(structure)
    rformula = get(info_dict, :formula, :None)
    num_fu = get(info_dict, :num_fu, 0)

    if (rformula != :None) & (num_fu != 0)
        return rformula, num_fu
    end

    # Find the reduced formula units
    sp_array = species(structure)
    unique_sp = sorted_symbols(unique(sp_array))
    num_atoms = Array{Int}(undef, size(unique_sp))
    for i in 1:length(unique_sp)
       num_atoms[i] = count(x -> x == unique_sp[i], sp_array)
    end
    num_fu = gcd(num_atoms)
    num_atoms ./= num_fu
    args = Symbol[]
    for i in 1:length(unique_sp)
        push!(args, unique_sp[i])
        if num_atoms[i] > 1
            push!(args, Symbol(num_atoms[i]))
        end
    end
    rformula = Symbol(args...)
    info(structure)[:formula] = rformula
    info(structure)[:num_fu] = num_fu
    return rformula, num_fu
end

"""
    specindex(structure::Cell)

Return the unique species and integer based indices for each Site 
"""
function specindex(structure)
    # Mapping between the speices as symbols and as intgers
    unique_spec = unique(species(structure))
    nunique = length(unique_spec)
    spec_indices = [findfirst(x-> x==sym, unique_spec) for sym in species(structure)]
    return unique_spec, spec_indices
end

function Base.show(io::IO, s::Cell)
    a, b, c, α, β, γ = cellpar(lattice(s))
    sym = join(map(string, species(structure)))
    print(io, "Cell $(sym) with $(nions(s)) atoms, lattice parametrs: $a $b $c $α $β $γ")
end

function Base.show(io::IO, ::MIME"text/plain", s::Cell)
    println(io, "Cell with $(nions(s)) ions")
    println(io, "Lattice: ")
    cellmat = s.lattice.matrix
    posmat = positions(s)
    for i in 1:3
        println(io, @sprintf "%8.3f   %8.3f   %8.3f"  cellmat[1, i] cellmat[2, i] cellmat[3, i])
    end

    println(io, "Sites: ")
    sym = species(s)
    for i in 1:nions(s)
        symbol = sym[i]
        println(io, @sprintf "%4s  %10.5f  %10.5f  %10.5f" symbol posmat[1, i] posmat[2, i] posmat[3, i])
    end
end

"""
    distance_matrix(structure::Cell; mic=true)

Compute the distance matrix for the given structure, using the minimum image convention (or not).
"""
function distance_matrix(structure::Cell)

    # Compute the naive pair-wise vectors
    nn = nions(structure)
    vecs = zeros(3, nn*nn)
    posmat = positions(structure)
    dmat = zeros(nn, nn)
    ivec = 0
    for i in 1:nn
        for j in i+1:nn
            ivec += 1
            for ii = 1:3
                vecs[ii, ivec] = posmat[ii, j] - posmat[ii, i]  
            end
        end
    end
    # Apply minimum image conventions
    _, dmic = mic(lattice(structure), @view(vecs[:, 1:ivec]))
    # Unpack computed distances
    ivec = 0
    for i in 1:nn
        for j in i+1:nn
            ivec += 1
            dmat[i, j] = dmic[ivec]
            dmat[j, i] = dmic[ivec]
        end
    end
    dmat
end


function distance_squared_between(posmat::Matrix, i, j, svec::Matrix, ishift)
    d2 = 0.
    for n=1:3
        d = posmat[n, j] - posmat[n, i] + svec[n, ishift]
        d2 += d *d
    end
    d2
end


"""
    check_minsep(structure::Strucuter, minsep::Dict)

Check if the minimum spearation constraints are satified. Minimum separations are supplied
as an dictionary where the global version under the :global key. To supply the minimum separations
between A and B, pass Dict((:A=>:B)=>1.0).
Return true or false.
"""
function check_minsep(structure::Cell, minsep::Dict{T, Float64}) where T
    _, spec_indices, minsep_matrix = compute_minsep_mat(structure, minsep)
    dist_mat = distance_matrix(structure)
    ni = nions(structure)
    for i in 1:ni
        for j in i+1:ni
            si = spec_indices[i]
            sj = spec_indices[j]
            if minsep_matrix[si, sj] > dist_mat[i, j]
                return false
            end
        end
    end
    return true
end


"""
    minsep_matrix(structure::Cell, minsep::Dict)

Initialis the minimum separation matrix and species mapping.
Returns the unique species, integer indexed sepcies and the minimum separation matrix.
"""
function compute_minsep_mat(structure, minsep::Dict{T, Float64}) where T
    unique_spec, spec_indices = specindex(structure)
    nunique = length(unique_spec)
    # Minimum separation matrix
    minsep_mat = zeros(nunique, nunique)
    for i in 1:nunique
        for j in i:nunique
            sA = unique_spec[i]
            sB = unique_spec[j]
            if (sA=>sB) in keys(minsep)
                value = minsep[sA=>sB]
            elseif (sB=>sA) in keys(minsep)
                value = minsep[sB=>sA]
            else
                value = get(minsep, :all, 1.0)
            end
            minsep_mat[i, j] = value
            minsep_mat[j, i] = value
        end
    end
    unique_spec, spec_indices, minsep_mat
end

"""
    make_supercell(structure::Cell, a, b, c)

Make a supercell

Currently only work with diagonal transform matrices.
TODO: Write function for the general cases....
"""
function make_supercell(structure::Cell, a, b, c)
    tmat = [
        a 0 0
        0 b 0
        0 0 c
    ]
    old_cell = cellmat(lattice(structure)) 
    new_cell = old_cell * tmat

    # Compute the shift vectors
    svec = shift_vectors(old_cell, a-1, b-1, c-1, 0, 0, 0)
    nshifts = size(svec)[2]

    current_pos = positions(structure)
    ns = natoms(structure)
    # New positions
    new_pos = zeros(3, nshifts * ns)
    for (i, shift) in enumerate(eachcol(svec))
        for j in 1:ns
            idx = j + (i - 1) * ns  # New index
            for n=1:3
                @inbounds new_pos[n, idx] = current_pos[n, j] + shift[n]
            end
        end
    end
    new_spec = repeat(species(structure), nshifts)
    Cell(Lattice(new_cell), new_pos, new_spec)
end


"""
    fingerprint(s::Cell; dmat=distance_matrix(s), weighted=true, cut_bl=3.0)

Computed the fingerprint vector based on simple sorted pair-wise distances.
"""
function fingerprint(s::Cell; dmat=distance_matrix(s), weighted=true, cut_bl=3.0)
    # Cut off distance based on minimum bond length
    cut_bl = minimum(d for d in dmat if d > 0.) * cut_bl   
    # Allocate workspace
    nn, _ = size(dmat)
    dist = zeros(nn * nn)

    # Weighting matrix
    num = atomic_numbers(s)

    c = 0
    norm_weights = 0.
    # Consider only the lower triangle
    for j in 1:nn
        for i in 1+j:nn
            d = dmat[i, j]
            d > cut_bl && continue
            c += 1
            dist[c] = d 
            if weighted
                dist[c] *= num[i] * num[j]
                norm_weights += num[i] * num[j]
            end
        end
    end

    dist_out = dist[1:c]

    # Normlise
    if weighted
        dist_out ./= (norm_weights / c)
    end
    sort!(dist_out)
end


"""
    fingerprint_distance(f1::AbstractVector, f2::AbstractVector;lim=Inf)

Compute the deviation between two finger print vectors

Comparison is truncated by the size of the shortest vector of the two, or by the `lim`
key word.
"""
function fingerprint_distance(f1::AbstractVector, f2::AbstractVector;lim=Inf)
    l1 = length(f1)
    l2 = length(f2)
    comp = min(l1, l2)
    d = 0.
    ncomp = 0
    for i in 1:comp
        f1[i] > lim && break 
        f2[i] > lim && break 
        d += abs(f1[i] - f2[i])
        ncomp += 1
    end
    d / comp
end