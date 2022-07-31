# Read dot castep files and exract trajectories 
module DotCastep

using LinearAlgebra
const NUM_TYPE = Float64

"""
Representation of a snapshot
"""
struct SnapShot{T}
    lattice::Matrix{T}
    species::Vector{Symbol}
    positions::Matrix{T}
    is_frac::Bool
    forces::Matrix{T}
    energy::T
    comment::Vector{String}
    extras::Dict
end

function SnapShot(lattice, species, positions, is_frac, force, energy, comment)
    SnapShot(lattice, species, positions, is_frac, force, energy, comment, Dict())
end



"""
    read_castep(fname::String)

Read a CASTEP file, return a Vector of the Snapshots
"""
function read_castep(fname::String; only_first=true)
    open(fname) do file
        read_castep(file; only_first=only_first, fname)
    end
end

function read_castep(io::IO; fname="Unkonwn file", only_first=true)
    lines = readlines(io)
    read_castep(lines; only_first, fname)
end



function read_castep(lines::Vector{String}; fname::AbstractString, only_first=true)
    nions = _get_nions(lines)
    nlines = length(lines)
    snapshots = SnapShot[]
    count = 1
    frame = 1
    while true
        if only_first
            count = skip_to_header(lines, count, nlines)
        end
        snapshot, count, ok, eof = _read_snapshot(lines, nions, count)
        eof && break
        !ok && continue
        snapshot.comment[] = "$fname - $frame"
        frame += 1
        push!(snapshots, snapshot)
    end
    return snapshots
end

"Skip to the next header"
function skip_to_header(lines::Vector{String}, count::Int, nlines::Int)
    if count > nlines
        return nlines
    end
    while true
        if occursin("CASTEP version", lines[count])
            break
        end
        count += 1
        if count > nlines
            break
        end
    end
    return count
end

"""
Get the number of ions
"""
function _get_nions(lines)
    nions = 0
    for line in lines
        if occursin("Total number of ions in cell", line)
            nions = parse(Int, split(line)[end])
            break
        end
    end
    return nions
end


"""
Read the initial structure and its energy/forces
"""
function _read_snapshot(lines::Vector{String}, nions::Int, offset::Int=0)
    count = 1 + offset

    forces = Array{NUM_TYPE}(undef, 3, nions)
    positions = Array{NUM_TYPE}(undef, 3, nions)
    stress = Array{NUM_TYPE}(undef, 3, 3)
    lattice = Array{NUM_TYPE}(undef, 3, 3)
    species = [:H for i in 1:nions]
    total_spin = 0.0

    coord_read = false
    eng_read = false
    force_read = false
    stress_read = false
    free_energy = 9999.9
    nlines = length(lines)
    eof = true
    while count < nlines
        line = lines[count]
        if occursin("Fractional coordinates of atoms", line)
            species = read_coord_table!(positions, @view lines[count+3:count+2+nions])
            coord_read = true
            count += nions + 2
        end
        if occursin("Integrated Spin Density", line)
            total_spin = parse(Float64, split(line)[5])
        end
        if occursin("Final free energy", line) & coord_read
            free_energy = parse(Float64, split(line)[end-1])
            eng_read = true
        end
        if occursin("Forces ********", line)
            read_coord_table!(forces, @view lines[count+6:count+5+nions])
            force_read = true
            count += nions + 5
            # Check if there are stress following the force to be read
            for subcount in count:(count+10)
                subline = lines[subcount]
                if occursin("Stress Tensor ********", subline)
                    read_coord_table!(stress, @view(lines[subcount+6:subcount+8]), 2, 1)
                    stress_read = true
                end
            end
            break
        end
        if occursin("Real Lattice(A)", line)
            lattice = read_lattice(@view lines[count+1:count+3])
        end
        count += 1
    end

    # Check if we have reached the end of file
    if count < nlines
        eof = false
    end

    abs_positions = lattice * positions
    extras = Dict{Symbol,Any}(:total_spin => total_spin,)
    if stress_read
        extras[:stress] = stress
    end
    snapshot = SnapShot(lattice, species, abs_positions, false, forces, free_energy, [""], extras)
    return snapshot, count, eng_read & force_read & coord_read, eof
end

"""
Read fractional coordinates (column vectors) from a table
"""
function read_coord_table!(property, lines, offset=3, specie_offset=2)
    species = Symbol[]
    for (iion, line) in enumerate(lines)
        tokens = split(line)
        if length(tokens) != (offset + 4)
            continue
        end
        property[1, iion] = parse(NUM_TYPE, tokens[offset+1])
        property[2, iion] = parse(NUM_TYPE, tokens[offset+2])
        property[3, iion] = parse(NUM_TYPE, tokens[offset+3])
        push!(species, Symbol(tokens[specie_offset]))
    end
    return species
end


"""
Read lattice vectors from the lines (column vectors)
"""
function read_lattice(lines)
    lattice = zeros(Float64, 3, 3)
    for i in 1:3
        tokens = split(lines[i])
        lattice[:, i] = map(x -> parse(NUM_TYPE, x), tokens[1:3])
    end
    return lattice
end

end # Module

using .DotCastep: read_castep

