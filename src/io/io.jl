#=
Code for input/output

Format:

- SHELX: read and write in AIRSS style
- CELL: read only
- XYZ: write only
- .castep: read only, including energies and forces
=#

include("io_cell.jl")
include("io_res.jl")
include("io_xyz.jl")
include("io_sheap.jl")
include("io_dotcastep.jl")

export read_res, write_res, read_cell, write_xyz
