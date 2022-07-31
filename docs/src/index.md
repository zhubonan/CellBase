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

```@autodocs
Modules = [CellBase]
Order = [:type, :function]
Pages = ["lattice.jl", "periodic.jl"]
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


## [Spglib.jl](https://github.com/singularitti/Spglib.jl) interface

The routine in `Spglib.jl` (which wraps the [spglib](https://spglib.github.io/spglib/) libraray) can be used for finding symmetry and performing reduction of the `Cell` type.

```@autodocs
Modules = [CellBase]
Pages = ["spg.jl"]
```


## Misc utils

Miscellaneous utility functions.

```@autodocs
Modules = [CellBase]
Order = [:type, :function]
Pages = ["minkowski.jl", "mathutils.jl"]
```


## Neighbour lists

Type and functions for handling neighbour lists.

```@autodocs
Modules = [CellBase]
Order = [:type, :function]
Pages = ["neighbour.jl"]
```

## IO

Code for read/writing file of crystal structures.

```@autodocs
Modules = [CellBase, CellBase.CellIO, CellBase.DotCastep, CellBase.SheapIO]
Order = [:type, :function]
Pages = ["io/io_cell.jl", "io/io_res.jl", "io/io_xyz.jl", "io/io_dotcastep.jl", "io/io_sheap.jl"]
```
