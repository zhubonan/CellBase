# CellBase.jl

Documentation for CellBase.jl

```@meta
CurrentModule = CellBase
```

# What is package does 

Provide a basic low level interface for storing and handling crystal structure data. 
The target application is for small and periodic cells, with the emphasis on both ease to use as well as performance.

# Guides

## `Lattice` type

The `Lattice` type represents the lattice.
The lattice vectors are stored as a matrix of **column vectors**. 
For example, the first lattice vector should be accessed as `cellmat(lattice)[:, 1]`.


!!! note "Updating an existing `Lattice`"

    The reciprocal lattice is also stored for quick access, however, this means that whenever the matrix is directly modified, the store reciprocal lattice vectors should be updated as well.
    Hence, the `set_cellmat!` function should be used for updating the cell matrix, which does this is automatically.


```@docs
Lattice
Lattice(matrix::Matrix{T}) where T
Lattice(va::Vector{T}, vb::Vector{T}, vc::Vector{T}) where T
cellmat
get_cellmat
reciprocal
rec_cellmat
get_rec_cellmat
update_rec!
cellvecs(lattice::Lattice)
volume(lattice::Lattice)
cellpar(lattice::Lattice)
wrap_positions(l::Lattice, v::AbstractVecOrMat)
mic
```

## `Cell` type

The `Cell` type represent a crystal structure combining the positions of atoms and the lattice vectors.
It basically combines the `Lattice`, the atomic positions and their identities.
Internally, the formere is stored as an matrix of column vectors of the absolute positions (cartesian space).

For flexibility, any additional array like data can be stored in the `arrays` property which is a dictionary with `Symbol` keys.
Additional metadata can be also be stored under `metadata`, which is a `Dict{Any, Any}`. 


```@autodocs
Modules = [CellBase]
Order = [:type, :function]
Pages = ["cell.jl"]
```