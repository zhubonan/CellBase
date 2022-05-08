#=
Interface for Spglib
=#
import Spglib
import Spglib:Cell as SCell

SCell(cell::Cell) = SCell(cellmat(cell), get_scaled_positions(cell), atomic_numbers(cell))

function Cell(cell::SCell)
    Cell(
        Lattice(collect(cell.lattice)),
        collect(cell.types),
        collect(cell.lattice * cell.positions),
    )
end

"Extend the Spglib methods"
macro extend_scell(func)
    quote
        Spglib.$(func)(cell::Cell, args...;kwargs...) = Spglib.$(func)(SCell(cell), args...;kwargs...)
    end
end

"Extend the Spglib methods that returns a Spglib.Cell"
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