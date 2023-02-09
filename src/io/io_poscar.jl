#=
Code for writing/parsing POSCAR
=#

"""
    write_poscar(io::IO, cell)

Write cell using the POSCAR format.
"""
function write_poscar(io::IO, cell::Cell)
    println(io, formula(cell))
    println(io, "  1.000000000000")
    latt = cellmat(cell)
    @printf io "    %20.16f %20.16f %20.16f\n" latt[1, 1] latt[2, 1] latt[3, 1]
    @printf io "    %20.16f %20.16f %20.16f\n" latt[1, 2] latt[2, 2] latt[3, 2]
    @printf io "    %20.16f %20.16f %20.16f\n" latt[1, 3] latt[2, 3] latt[3, 3]

    # Write species
    # Here we need to sort the cell
    syms = species(cell)
    usyms = unique(syms)
    # Array of symbol index
    sindex = [findfirst(x -> x == sym, usyms) for sym in syms]

    sort_idx = sortperm(sindex)
    println(io, "   ", join(map(string, usyms), "  "))

    # Check the number of species
    specie_counts = [count(x -> x == sym, syms) for sym in usyms]
    println(io, "   ", join(map(string, specie_counts), "  "))
    println(io, "Direct")
    scaled_pos = get_scaled_positions(cell)
    for idx in sort_idx
        x, y, z = scaled_pos[:, idx]
        @printf io "  %20.16f  %20.16f  %20.16f\n" x y z
    end

end

"""
    write_poscar(fname::AbstractString, cell::Cell)
Write `cell` as a `POSCAR` file.
"""
function write_poscar(fname::AbstractString, cell::Cell)
    open(fname, "w") do io
        write_poscar(io, cell)
    end
end

export write_poscar
