#=
Functions for reading SHELX files
=#

using Printf
"""
Read an array containing the lines of the SHELX file
"""
function read_res(lines::Vector{String})
    cellpar = Array{Float64}(undef, 6)
    title_items = Dict{Symbol, Any}()
    line_no = 1 

    species = fill(:na, length(lines))
    scaled_pos = zeros(3, length(lines))
    spins = zeros(3, length(lines))
    iatom = 1
    while line_no < length(lines)
        line = lines[line_no]
        tokens = split(strip(line))
        if tokens[1] == "TITL" 
            title_items = parse_titl(line)
        elseif (tokens[1] == "CELL") & (length(tokens) == 8)
            cellpar[:] = map(x -> parse(Float64, x), tokens[3:8])
        elseif tokens[1] == "SFAC"
            for atom_line in @view lines[line_no+1:end]
                if strip(atom_line) == "END"
                    break
                end
                atokens = split(strip(atom_line))
                species[iatom] = Symbol(atokens[1])
                scaled_pos[:, iatom] = parse.(Float64, atokens[3:5])
                if length(atokens) == 7
                    spins[iatom] = parse(Float64, atokens[7])
                else
                    spins[iatom] = 0.
                end
                iatom += 1
            end
        end
        line_no += 1
    end
    # Adjust the sizes
    spins = spins[1:iatom-1]
    scaled_pos = @view scaled_pos[:, 1:iatom-1]
    species = species[1:iatom-1]

    lattice = Lattice(cellpar)
    cell = Cell(lattice, species, cellmat(lattice) * scaled_pos)

    # Attach spin only if there are any non-zero ones...
    if any(x ->x != 0, spins) 
        cell.arrays[:spins] = spins
    end

    CellBase.attachmetadata!(cell, title_items)
    cell
end

function read_res(s::AbstractString)
    open(s) do handle
        read_res(readlines(handle))
    end
end

function parse_titl(s::AbstractString)
    pfloat(x) = parse(Float64, x)
    tokens = split(strip(s))[2:end]
    Dict(
        :label=>tokens[1],
        :pressure=>pfloat(tokens[2]),
        :volume=>pfloat(tokens[3]),
        :enthalpy=>pfloat(tokens[4]),
        :spin=>pfloat(tokens[5]),
        :abs_spin=>pfloat(tokens[6]),
        :natoms=>pfloat(tokens[7]),
        :symm=>tokens[8],
        :flag1=> tokens[9],
        :flag2=> tokens[10],
        :flag3=> tokens[11],
    )
end

"""
    write_res(io::IO, structure::Cell)

Write out SHELX format data
"""
function write_res(io::IO, structure::Cell)
    infodict = structure.metadata
    titl = (
        label=get(infodict, :label, "jurss-in-out"),
        pressure=get(infodict, :pressure, 0.0),
        volume=volume(structure),
        enthalpy=get(infodict, :enthalpy, 0.0),
        spin=get(infodict, :spin, 0.0),
        abs_spin=get(infodict, :abs_spin, 0.0),
        natoms=nions(structure),
        symm=get(infodict, :symm, "(n/a)"),
        flag1= get(infodict, :flag1, "n"),
        flag2= get(infodict, :flag2, "-"),
        flag3= get(infodict, :flag3, "1")
    ) 
    titl_line = @sprintf("TITL %s %.10f %.10f %.10f %.3f %.3f %d %s %s %s %s\n", titl...)
    write(io, titl_line)
    cell_line = @sprintf("CELL 1.54180 %.6f %.6f %.6f %.6f %.6f %.6f\n", cellpar(lattice(structure))...)
    write(io, cell_line)
    write(io, "LATT -1\n")
    write(io, "SFAC ", join(map(string, unique(species(structure))), " "), "\n")
    fposmat = CellBase.get_scaled_positions(structure)
    # Wrap
    fposmat .-= floor.(fposmat)

    count = 1
    last_symbol = structure.symbols[1]
    if :spins in keys(structure.arrays)
        spin_array = structure.arrays[:spins]
        for (i, symbol) in enumerate(structure.symbols)
            if symbol != last_symbol
                count +=1
            end
            write(io, @sprintf("%-7s %2s %15.12f %15.12f %15.12f 1.0 %.3f\n", symbol, count, fposmat[1, i], fposmat[2, i], fposmat[3, i], spin_array[i]))
            last_symbol = symbol
        end
    else
        for (i, symbol) in enumerate(structure.symbols)
            if symbol != last_symbol
                count +=1
            end
            write(io, @sprintf("%-7s %2s %15.12f %15.12f %15.12f 1.0\n", symbol, count, fposmat[1, i], fposmat[2, i], fposmat[3, i]))
            last_symbol = symbol
        end
    end
    write(io, "END\n")
end

"""
    write_res(fname::AbstractString, structure::Cell)

Write out SHELX format data to a file
"""
function write_res(fname::AbstractString, structure::Cell)
    open(fname, "w") do fh 
        write_res(fh, structure)
    end
end