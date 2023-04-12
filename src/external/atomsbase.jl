import AtomsBase
using Unitful


function _cell_to_atomic_system(
    lattice::Lattice{T},
    species::AbstractArray,
    positions::AbstractArray,
    arrays::Dict;
    cell_unit=u"Å",
    kwargs...,
) where {T}
    cm = cellmat(lattice)
    ndim = size(cm, 1)
    box = SVector{ndim,SVector{ndim,T}}(SVector{ndim,T}(col) for col in eachcol(cm))
    # Unpack :arrays
    kws = [Dict{Symbol,Any}() for _ = 1:length(species)]
    for (key, value) in arrays
        for (iat, prop) in enumerate(eachslice(value; dims=ndims(value)))
            kws[iat][key] = prop
        end
    end

    atoms = [
        AtomsBase.Atom(symbol, pos .* cell_unit; kw...) for
        (symbol, pos, kw) in zip(species, eachcol(positions), kws)
    ]
    AtomsBase.periodic_system(atoms, box .* cell_unit; kwargs...)
end

"""
    atomic_system(cell::Cell;cell_unit=u"Å")

Construct a `Cell` object into a `FlexibleSystem`.
The length unit for the `Cell` object passed is assuemd to be `cell_unit`.
"""
AtomsBase.atomic_system(cell::Cell; cell_unit=u"Å") = _cell_to_atomic_system(
    lattice(cell),
    species(cell),
    positions(cell),
    cell.arrays;
    cell_unit,
    metadata(cell)...,
)

"""
    Cell(system::AbstractSystem;cell_unit=u"Å")

Construct a `Cell` object from AbstractSystem.
The length unit for the `Cell` object returned well be in the `cell_unit`.
"""
function Cell(system::AtomsBase.AbstractSystem; cell_unit=u"Å")
    pos = hcat(map(x -> collect(ustrip.(cell_unit, x)), AtomsBase.position(system))...)
    cm = hcat(map(x -> collect(ustrip.(cell_unit, x)), AtomsBase.bounding_box(system))...)
    @assert all(AtomsBase.periodicity(system))
    out = Cell(Lattice(cm), AtomsBase.atomic_symbol(system), pos)
    # Store sys.data in metadata
    if isa(system, AtomsBase.FlexibleSystem)
        for (key, value) in system.data
            out.metadata[key] = value
        end
        property_keys = vcat(map(x -> collect(keys(x.data)), system)...) |> unique
        for key in property_keys
            out.arrays[key] = hcat(map(x -> x[key], system)...)
        end
    end
    out
end
