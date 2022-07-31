#=
Functions for writing XYZ files
=#

## XYZ files

"""
    Write snapshots to a xyz file
"""
function write_xyz(fname, structures::Vector{Cell{T}}) where T
    lines = String[]
    for structure in structures
        push_xyz!(lines, structure)
    end
    
    open(fname, "w") do handle
        for line in lines
            write(handle, line)
            write(handle, "\n")
        end
    end
end


"""xyz lines for a single frame"""
function push_xyz!(lines, structure::Cell)
    ns = nions(structure)
    cell = cellmat(lattice(structure))
    push!(lines, "$ns")
    ax, ay, az = cell[:, 1]
    bx, by, bz = cell[:, 2]
    cx, cy, cz = cell[:, 3]

    # Include extra properties
    info_lines = []
    for (key, value) in CellBase.metadata(structure)
        push!(info_lines, "$(key)=\"$(value)\"")
    end
    info_string = join(info_lines, " ")
    
    comment_line = "Lattice = \"$ax $ay $az $bx $by $bz $cx $cy $cz\" Properties=\"species:S:1:pos:R:3\" $(info_string)"
    push!(lines, comment_line)
    
    # Write the atoms
    sp = species(structure)
    pos = positions(structure)
    for i in 1:ns
        specie = sp[i]
        x, y, z = pos[:, i]
        line = "$specie $x $y $z"
        push!(lines, line)
    end
    return lines
end