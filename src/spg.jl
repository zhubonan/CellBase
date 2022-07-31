#=
Interface for Spglib
=#


import Spglib
import Spglib: get_symmetry, get_international, get_dataset

"Alias for `Spglib.Cell`"
const SCell = Spglib.Cell

"""
    SCell(cell::Cell)

Construct `Spglib.Cell` from `Cell` type.
"""
SCell(cell::Cell) = SCell(cellmat(cell), get_scaled_positions(cell), atomic_numbers(cell))

"""
    Cell(cell::SCell)

Return a `Cell` object from `Spglib.Cell`.
"""
function Cell(cell::SCell)
    Cell(
        Lattice(collect(cell.lattice)),
        collect(cell.types),
        collect(cell.lattice * cell.positions),
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
        Spglib.$(func)(cell::Cell, args...;kwargs...) = Cell(Spglib.$(func)(SCell(cell), args...;kwargs...))
    end
end

# Extend 
@extend_scell get_dataset
@extend_scell get_symmetry
@extend_scell get_international
@extend_scell_roundtrip standardize_cell
@extend_scell_roundtrip find_primitive
@extend_scell_roundtrip refine_cell
@extend_scell_roundtrip niggli_reduce
@extend_scell_roundtrip delaunay_reduce

export get_symmetry, get_dataset, get_international, refine_cell, niggli_reduce, delaunay_reduce