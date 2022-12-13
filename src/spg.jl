import Base
import Spglib
import Spglib: get_symmetry, get_international, get_dataset
import Spglib: standardize_cell, find_primitive, refine_cell, SpglibError

struct SpglibConversionError <: Exception
    var::String
end
Base.showerror(io::IO, e::SpglibConversionError) = print(io, "Cannot convert to/from with Spglib.jl types consistently:", e.var, "!")

"Alias for `Spglib.Cell`"
const SCell = Spglib.Cell

"""
    SCell(cell::Cell)

Construct `Spglib.Cell` from `Cell` type.
"""
SCell(cell::Cell) = SCell(cellmat(cell), collect(eachcol(get_scaled_positions(cell))),
                          atomic_numbers(cell))

"""
    Cell(cell::SCell)

Return a `Cell` object from `Spglib.Cell`.
"""
function Cell(cell::SCell)
    # Collect positions
    pos = zeros(Float64, 3, length(cell.positions))
    for i in 1:length(cell.positions)
        pos[:, i] .= cell.lattice * cell.positions[i]
    end
    Cell(
        Lattice(collect(cell.lattice)),
        collect(cell.types),
        pos
    )
end

"""
Macro for extending the Spglib methods.

Usage:

```
@extend_scell get_dataset
```

will allow the `get_dataset` method of Spglib to be used for `Cell` type.
"""
macro extend_scell(func)
    quote
        Spglib.$(func)(cell::Cell, args...;kwargs...) = Spglib.$(func)(SCell(cell), args...;kwargs...)
    end
end

"""
Macro for extending the Spglib methods and convert returned `Spglib.Cell` to `Cell`.

Usage:

```
@extend_scell_roundtrip standardize_cell
```

will allow the `standardize_cell` method to be used and the returned `Spglib.Cell` is converted to `Cell`.
"""
macro extend_scell_roundtrip(func)
    quote
        function Spglib.$(func)(cell::Cell, args...;kwargs...)
            new_cell = Cell(Spglib.$(func)(SCell(cell), args...;kwargs...))

            # Create a copy of the original cell and set positons and cell
            out_cell = deepcopy(cell)
            set_cellmat!(out_cell, cellmat(new_cell))
            set_positions!(out_cell, positions(new_cell))
            # Check for any reordering of the symbols
            if any(new_cell.symbols .!= out_cell.symbols)
                if !isempty(out_cell.arrays)
                    throw(SpglibConversionError("atomic ordering has changed!"))
                else
                    @warn "Order of the symbols has changed!"
                    out_cell.symbols = new_cell.symbols
                end
            end
            out_cell
        end
    end
end

# Extend 
@extend_scell get_dataset
@extend_scell get_symmetry
@extend_scell get_international
@extend_scell get_multiplicity
@extend_scell get_dataset_with_hall_number
@extend_scell get_spacegroup_number
@extend_scell get_spacegroup_type
@extend_scell get_schoenflies

@extend_scell_roundtrip standardize_cell
@extend_scell_roundtrip find_primitive
@extend_scell_roundtrip refine_cell
@extend_scell_roundtrip niggli_reduce
@extend_scell_roundtrip delaunay_reduce

"""
    niggli_reduce_cell(cell::Cell; wrap_pos=true)

Apply niggli reduction to the lattice using Spglib 
"""
function niggli_reduce_cell(cell::Cell; wrap_pos=true)
    reduced_cellmat = Spglib.niggli_reduce(cellmat(cell))
    outcell = deepcopy(cell)
    set_cellmat!(outcell, reduced_cellmat, scale_positions=false)
    wrap_pos && wrap!(outcell)
    outcell
end

export get_symmetry, get_dataset, get_international, refine_cell, standardize_cell, find_primitive, niggli_reduce_cell