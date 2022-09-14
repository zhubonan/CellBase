var documenterSearchIndex = {"docs":
[{"location":"#CellBase.jl","page":"CellBase.jl","title":"CellBase.jl","text":"","category":"section"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"Documentation for CellBase.jl","category":"page"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"CurrentModule = CellBase","category":"page"},{"location":"#What-is-package-does","page":"CellBase.jl","title":"What is package does","text":"","category":"section"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"Provide a basic low level interface for storing and handling crystal structure data.  The target application is for small and periodic cells, with the emphasis on both ease to use as well as performance.","category":"page"},{"location":"#Guides","page":"CellBase.jl","title":"Guides","text":"","category":"section"},{"location":"#Lattice-type","page":"CellBase.jl","title":"Lattice type","text":"","category":"section"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"The Lattice type represents the lattice. The lattice vectors are stored as a matrix of column vectors.  For example, the first lattice vector should be accessed as cellmat(lattice)[:, 1].","category":"page"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"note: Updating an existing `Lattice`\nThe reciprocal lattice is also stored for quick access, however, this means that whenever the matrix is directly modified, the store reciprocal lattice vectors should be updated as well. Hence, the set_cellmat! function should be used for updating the cell matrix, which does this is automatically.","category":"page"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"Modules = [CellBase]\nOrder = [:type, :function]\nPages = [\"lattice.jl\", \"periodic.jl\"]","category":"page"},{"location":"#CellBase.Lattice","page":"CellBase.jl","title":"CellBase.Lattice","text":"Lattice{T}\n\nThe type represents a lattice in 3D space formed by three column lattice vectors.\n\n\n\n\n\n","category":"type"},{"location":"#CellBase.Lattice-Union{Tuple{Matrix{T}}, Tuple{T}} where T","page":"CellBase.jl","title":"CellBase.Lattice","text":"Lattice(matrix::Matrix{T}) where T\n\nConstruct a Lattice from a matrix of column vectors.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.Lattice-Union{Tuple{T}, NTuple{6, T}} where T","page":"CellBase.jl","title":"CellBase.Lattice","text":"Lattice(a::T, b::T, c::T, α::T, β::T, γ::T) where T\n\nConstruct a Lattice from lattice parameters.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.Lattice-Union{Tuple{T}, Tuple{T, T, T}} where T<:Real","page":"CellBase.jl","title":"CellBase.Lattice","text":"Lattice(a::T, b::T, c::T) where T <: Real\n\nConstruct a Lattice with orthogonal lattice vectors.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.Lattice-Union{Tuple{T}, Tuple{Vector{T}, Vector{T}, Vector{T}}} where T","page":"CellBase.jl","title":"CellBase.Lattice","text":"Lattice(va::Vector{T}, vb::Vector{T}, vc::Vector{T}) where T\n\nConstruct a Lattice from three lattice vectors.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.Lattice-Union{Tuple{Vector{T}}, Tuple{T}} where T","page":"CellBase.jl","title":"CellBase.Lattice","text":"Lattice(cellpar::Vector{T}) where T\n\nConstruct a Lattice from a six-vector of lattice parameters.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.cellmat-Tuple{Lattice}","page":"CellBase.jl","title":"CellBase.cellmat","text":"cellmat(lattice::Lattice)\n\nReturn the matrix of column vectors.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.cellmat_row-Tuple{Lattice}","page":"CellBase.jl","title":"CellBase.cellmat_row","text":"Get a matrix of column vectors\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.cellpar-Tuple{Lattice}","page":"CellBase.jl","title":"CellBase.cellpar","text":"cellpar(lattice::Lattice)\n\nReturn the lattice parameters as a six-vector.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.cellvecs-Tuple{Lattice}","page":"CellBase.jl","title":"CellBase.cellvecs","text":"Lattice vectors\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.frac_pos-Tuple{Site, Lattice}","page":"CellBase.jl","title":"CellBase.frac_pos","text":"Fraction positions of a site\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.get_cellmat-Tuple{Lattice}","page":"CellBase.jl","title":"CellBase.get_cellmat","text":"Return a copy of the cell matrix\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.get_rec_cellmat-Tuple{Lattice}","page":"CellBase.jl","title":"CellBase.get_rec_cellmat","text":"Return a copy of the reciprocal cell matrix\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.get_scaled_positions-Tuple{Lattice, AbstractVecOrMat{T} where T}","page":"CellBase.jl","title":"CellBase.get_scaled_positions","text":"get_scaled_positions(l::Lattice, v::AbstractVecOrMat) = reciprocal(l) * v\n\nReturn the scaled positions.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.mic-Tuple{Lattice, AbstractVecOrMat{T} where T}","page":"CellBase.jl","title":"CellBase.mic","text":"mic(l::Lattice, vec::AbstractVecOrMat)\n\nCompute the minimum-image convention representation of a series of displacement vectors.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.mic_naive-Tuple{Lattice, AbstractMatrix{T} where T}","page":"CellBase.jl","title":"CellBase.mic_naive","text":"Compute naive MIC representation of the vector(s)\n\nNot safe for skewed cells. Requires\n\nnorm(length) < 0.5 * min(cellpar(l)[1:3])\n\nSee also:\n\nW. Smith, \"The Minimum Image Convention in Non-Cubic MD Cells\", 1989,   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.57.1696.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.mic_naive-Tuple{Lattice, AbstractVector{T} where T}","page":"CellBase.jl","title":"CellBase.mic_naive","text":" mic_safe(l::Lattice, v::AbstractVector)\n\nCompute MIC representation of vectors using a safe approach based on Minkowski reduced cell\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.mic_safe-Tuple{Lattice, AbstractMatrix{T} where T}","page":"CellBase.jl","title":"CellBase.mic_safe","text":" mic_safe(l::Lattice, v::AbstractMatrix)\n\nCompute MIC representation of vectors using a safe approach based on Minkowski reduced cell\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.mic_safe-Tuple{Lattice, AbstractVector{T} where T}","page":"CellBase.jl","title":"CellBase.mic_safe","text":" mic_safe(l::Lattice, v::AbstractVector)\n\nCompute MIC representation of vectors using a safe approach based on Minkowski reduced cell\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.random_vec-Tuple{Lattice}","page":"CellBase.jl","title":"CellBase.random_vec","text":"Random displacement vector in a unit cell\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.random_vec_in_cell-Union{Tuple{Matrix{T}}, Tuple{T}} where T","page":"CellBase.jl","title":"CellBase.random_vec_in_cell","text":"random_vec_in_cell(cell::Matrix{T}; scales::Vector=ones(size(cell)[1])) where T\n\nGet an random vector within a unit cell. The cell is a matrix made of column vectors.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.rec_cellmat-Tuple{Lattice}","page":"CellBase.jl","title":"CellBase.rec_cellmat","text":"rec_cellmat(l::Lattice)\n\nReturns the matrix of reciprocal lattice vectors (without the 2pi factor).\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.reciprocal-Tuple{Lattice}","page":"CellBase.jl","title":"CellBase.reciprocal","text":"reciprocal(l::Lattice)\n\nReturns the matrix of reciprocal lattice vectors (without the 2pi factor).\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.scaled_positions","page":"CellBase.jl","title":"CellBase.scaled_positions","text":"Alias  to get_scaled_positions\n\n\n\n\n\n","category":"function"},{"location":"#CellBase.set_cellmat!-Tuple{Lattice, Any}","page":"CellBase.jl","title":"CellBase.set_cellmat!","text":"set_cellmat!(lattice::Lattice, mat)\n\nSet the cell matrix and update the reciprocal cell matrix.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.update_rec!-Tuple{Lattice}","page":"CellBase.jl","title":"CellBase.update_rec!","text":"update_rec!(l::Lattice)\n\nUpdate the reciprocal matrix - this should be called everytime cell is changed.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.volume-Tuple{Lattice}","page":"CellBase.jl","title":"CellBase.volume","text":"Return the volume of the cell\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.wrap!-Tuple{Site, Lattice}","page":"CellBase.jl","title":"CellBase.wrap!","text":"Wrap a site back into the box\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.wrap_positions-Tuple{Lattice, AbstractVecOrMat{T} where T}","page":"CellBase.jl","title":"CellBase.wrap_positions","text":"wrap_positions(l::Lattice, v::AbstractVecOrMat)\n\nWrap positions back into the box defined by the cell vectors.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.max_shifts-Tuple{AbstractArray, Real}","page":"CellBase.jl","title":"CellBase.max_shifts","text":"Compute the require periodicities required to fill a cut off sphere of a certain radius\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.shift_indices-Tuple{AbstractArray, Real}","page":"CellBase.jl","title":"CellBase.shift_indices","text":"shift_indices(lattice::AbstractArray, rc::Real)\n\nCompute a vector containing shift indices\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.shift_vectors","page":"CellBase.jl","title":"CellBase.shift_vectors","text":"shift_vectors(lattice::AbstractMatrix, shift1, shift2, shift3)\n\nCompute the shift vectors needed to include all parts of the lattice within a cut off radius.\n\n\n\n\n\n","category":"function"},{"location":"#CellBase.shift_vectors-Tuple{AbstractMatrix{T} where T, Real}","page":"CellBase.jl","title":"CellBase.shift_vectors","text":"Compute shift vectors given lattice matrix and cut off radius\n\n\n\n\n\n","category":"method"},{"location":"#Cell-type","page":"CellBase.jl","title":"Cell type","text":"","category":"section"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"The Cell type represent a crystal structure combining the positions of atoms and the lattice vectors. It basically combines the Lattice, the atomic positions and their identities. Internally, the formere is stored as an matrix of column vectors of the absolute positions (cartesian space).","category":"page"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"For flexibility, any additional array like data can be stored in the arrays property which is a dictionary with Symbol keys. Additional metadata can be also be stored under metadata, which is a Dict{Any, Any}. ","category":"page"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"Modules = [CellBase]\nOrder = [:type, :function]\nPages = [\"cell.jl\"]","category":"page"},{"location":"#CellBase.Cell","page":"CellBase.jl","title":"CellBase.Cell","text":"A Cell represents a periodic structure in three-dimensional space.\n\nDefined as:\n\nmutable struct Cell{T}\n    lattice::Lattice{T}                 # Lattice of the structure\n    symbols::Vector{Symbol}\n    positions::Matrix{T}\n    arrays::Dict{Symbol, AbstractArray}        # Any additional arrays\n    metadata::Dict\nend\n\n\n\n\n\n","category":"type"},{"location":"#CellBase.Cell-Tuple{Lattice, Any, Vector{T} where T}","page":"CellBase.jl","title":"CellBase.Cell","text":"Cell(lat::Lattice, numbers::Vector{Int}, positions::Vector)\n\nConstructure the Cell type from lattice, positions and numbers\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.Cell-Tuple{Lattice, Vector{Symbol}, Matrix{T} where T}","page":"CellBase.jl","title":"CellBase.Cell","text":"Cell(l::Lattice, symbols, positions) where T\n\nConstruct a Cell type from arrays\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.Cell-Union{Tuple{T}, Tuple{Lattice, Vector{T}, Matrix{T} where T}} where T<:Real","page":"CellBase.jl","title":"CellBase.Cell","text":"Cell(lat::Lattice, numbers::Vector{Int}, positions)\n\nConstructure the Cell type from lattice, positions and numbers\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.array-Tuple{Cell, Symbol}","page":"CellBase.jl","title":"CellBase.array","text":"array(structure::Cell, arrayname::Symbol)\n\nReturn the additional array stored in the Cell object.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.arraynames-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.arraynames","text":"arraynames(structure::Cell)\n\nReturn the names of additional arrays.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.atomic_numbers-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.atomic_numbers","text":"atomic_numbers(structure::Cell)\n\nReturn a Vector of the atomic numbers.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.attachmetadata!-Tuple{Cell, Dict}","page":"CellBase.jl","title":"CellBase.attachmetadata!","text":"attachmetadata!(structure::Cell, metadata::Dict)\n\nReplace metadata with an existing dictionary.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.cellmat-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.cellmat","text":"cellmat(structure::Cell)\n\nReturn the matrix of lattice vectors.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.cellpar-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.cellpar","text":"cellpar(structure::Cell)\n\nReturn the lattice parameters.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.check_minsep-Union{Tuple{T}, Tuple{Cell, Dict{T, Float64}}} where T","page":"CellBase.jl","title":"CellBase.check_minsep","text":"check_minsep(structure::Cell, minsep::Dict)\n\nCheck if the minimum separation constraints are satisfied. Minimum separations are supplied as an dictionary where the global version under the :global key. To supply the minimum separations between A and B, pass Dict((:A=>:B)=>1.0). Return true or false.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.clip-Tuple{Cell, AbstractVector{T} where T}","page":"CellBase.jl","title":"CellBase.clip","text":"clip(s::Cell, mask::AbstractVector)\n\nClip a structure with a given indexing array\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.compute_minsep_mat-Union{Tuple{T}, Tuple{Any, Dict{T, Float64}}} where T","page":"CellBase.jl","title":"CellBase.compute_minsep_mat","text":"minsep_matrix(structure::Cell, minsep::Dict)\n\nInitialise the minimum separation matrix and species mapping. Returns the unique species, integer indexed species and the minimum separation matrix.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.distance_matrix-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.distance_matrix","text":"distance_matrix(structure::Cell; mic=true)\n\nCompute the distance matrix for the given structure, using the minimum image convention (or not).\n\nNote the returned matrix does not expand the cell. The distance matrix cannot be safety used for obtaining the minimum separations. For example, a structure with a single atom would be a distance matrix containing only zero.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.distance_squared_between-Tuple{Matrix{T} where T, Any, Any, Matrix{T} where T, Any}","page":"CellBase.jl","title":"CellBase.distance_squared_between","text":"distance_squared_between(posmat::Matrix, i, j, svec::Matrix, ishift)\n\nReturn the squared distance between two positions stored in a matrix and shift vector.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.fingerprint-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.fingerprint","text":"fingerprint(s::Cell; dmat=distance_matrix(s), weighted=true, cut_bl=3.0)\n\nComputed the fingerprint vector based on simple sorted pair-wise distances. NOTE: Does not work for single atom cell!!\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.fingerprint_distance-Tuple{AbstractVector{T} where T, AbstractVector{T} where T}","page":"CellBase.jl","title":"CellBase.fingerprint_distance","text":"fingerprint_distance(f1::AbstractVector, f2::AbstractVector;lim=Inf)\n\nCompute the deviation between two finger print vectors\n\nComparison is truncated by the size of the shortest vector of the two, or by the lim key word.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.formula_and_factor-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.formula_and_factor","text":"formula_and_factor(structure::Cell)\n\nReturn reduced formula and associated factor.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.get_cellmat-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.get_cellmat","text":"get_cellmat(structure::Cell)\n\nReturn the matrix of lattice vectors(copy).\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.get_lattice-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.get_lattice","text":"get_lattice(structure::Cell)\n\nReturn the Lattice instance (copy).\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.get_positions-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.get_positions","text":"get_positions(cell::Cell)\n\nReturn a copy of the positions (cartesian coordinates) of the atoms in a structure.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.get_scaled_positions-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.get_scaled_positions","text":"get_fraction_positions(cell::Cell)\n\nReturn fractional positions of the cell.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.get_species-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.get_species","text":"get_species(structure::Cell)\n\nReturn a Vector (copy) of species names.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.lattice-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.lattice","text":"lattice(structure::Cell)\n\nReturn the Lattice instance.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.make_supercell-Tuple{Cell, Any, Any, Any}","page":"CellBase.jl","title":"CellBase.make_supercell","text":"make_supercell(structure::Cell, a, b, c)\n\nMake a supercell\n\nCurrently only work with diagonal transform matrices. TODO: Write function for the general cases....\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.metadata-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.metadata","text":"metadata(structure::Cell)\n\nReturn the metadata dictionary. \n\n\n\n\n\n","category":"method"},{"location":"#CellBase.natoms","page":"CellBase.jl","title":"CellBase.natoms","text":"natoms(cell::Cell)\n\nReturn number of atoms in a structure.\n\n\n\n\n\n","category":"function"},{"location":"#CellBase.nions-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.nions","text":"nions(cell::Cell)\n\nReturn number of atoms in a structure.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.num_fu-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.num_fu","text":"num_fu(structure::Cell)\n\nReturn the number of formula units.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.positions-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.positions","text":"positions(cell::Cell)\n\nReturn positions (cartesian coordinates) of the atoms in a structure.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.rattle!-Tuple{Cell, Any}","page":"CellBase.jl","title":"CellBase.rattle!","text":"rattle!(cell::Cell, amp)\n\nRattle the positions of the cell for a given maximum amplitude (uniform distribution).\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.reduced_fu-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.reduced_fu","text":"reduced_fu(structure::Cell)\n\nReturn the reduced formula.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.set_cellmat!-Tuple{Cell, Any}","page":"CellBase.jl","title":"CellBase.set_cellmat!","text":"set_cellmat!(cell::Cell, mat;scale_positions=true)\n\nUpdate the Lattice with a new matrix of lattice vectors.  Scale of the existing postions if needed.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.set_positions!-Tuple{Cell, Any}","page":"CellBase.jl","title":"CellBase.set_positions!","text":"set_positions!(cell::Cell, pos)\n\nSet the positions of the Cell with a new matrix.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.set_scaled_positions!-Tuple{Cell, Matrix{T} where T}","page":"CellBase.jl","title":"CellBase.set_scaled_positions!","text":"set_scaled_positions!(cell::Cell, scaled::Matrix)\n\nSet scaled positions for a cell.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.sorted_symbols-Tuple{Any}","page":"CellBase.jl","title":"CellBase.sorted_symbols","text":"sorted_symbols(symbols)\n\nReturn sorted symbols by atomic numbers.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.species-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.species","text":"species(structure::Cell)\n\nReturn a Vector of species names.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.specindex-Tuple{Any}","page":"CellBase.jl","title":"CellBase.specindex","text":"specindex(structure::Cell)\n\nReturn the unique species and integer based indices for each atom.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.sposarray-Union{Tuple{Cell{T}}, Tuple{T}} where T","page":"CellBase.jl","title":"CellBase.sposarray","text":"sposarray(cell::Cell)\n\nReturn the positions as a Vector of static arrays. The returned array can provide improved performance for certain type of operations.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.volume-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.volume","text":"volume(structure::Cell)\n\nReturn the volume of the cell.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.wrap!-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.wrap!","text":"wrap!(cell::Cell)\n\nWrap atom outside of the lattice back into the box defined by the lattice vectors.\n\n\n\n\n\n","category":"method"},{"location":"#[Spglib.jl](https://github.com/singularitti/Spglib.jl)-interface","page":"CellBase.jl","title":"Spglib.jl interface","text":"","category":"section"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"The routine in Spglib.jl (which wraps the spglib libraray) can be used for finding symmetry and performing reduction of the Cell type.","category":"page"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"Modules = [CellBase]\nPages = [\"spg.jl\"]","category":"page"},{"location":"#CellBase.Cell-Tuple{Spglib.Cell}","page":"CellBase.jl","title":"CellBase.Cell","text":"Cell(cell::SCell)\n\nReturn a Cell object from Spglib.Cell.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.SCell","page":"CellBase.jl","title":"CellBase.SCell","text":"Alias for Spglib.Cell\n\n\n\n\n\n","category":"type"},{"location":"#CellBase.SCell-Tuple{Cell}","page":"CellBase.jl","title":"CellBase.SCell","text":"SCell(cell::Cell)\n\nConstruct Spglib.Cell from Cell type.\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.@extend_scell-Tuple{Any}","page":"CellBase.jl","title":"CellBase.@extend_scell","text":"Macro for extending the Spglib methods.\n\nUsage:\n\n@extend_scell get_dataset\n\nwill allow the get_dataset method of Spglib to be used for Cell type.\n\n\n\n\n\n","category":"macro"},{"location":"#CellBase.@extend_scell_roundtrip-Tuple{Any}","page":"CellBase.jl","title":"CellBase.@extend_scell_roundtrip","text":"Macro for extending the Spglib methods and convert returned Spglib.Cell to Cell.\n\nUsage:\n\n@extend_scell_roundtrip standardize_cell\n\nwill allow the standardize_cell method to be used and the returned Spglib.Cell is converted to Cell.\n\n\n\n\n\n","category":"macro"},{"location":"#Misc-utils","page":"CellBase.jl","title":"Misc utils","text":"","category":"section"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"Miscellaneous utility functions.","category":"page"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"Modules = [CellBase]\nOrder = [:type, :function]\nPages = [\"minkowski.jl\", \"mathutils.jl\"]","category":"page"},{"location":"#CellBase.reduction_gauss-Tuple{Matrix{T} where T, Vector{T} where T, Vector{T} where T}","page":"CellBase.jl","title":"CellBase.reduction_gauss","text":"2D Gauss reduction\n\nWorks for column vectors instead of row vectors\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.cellpar2mat-NTuple{6, Any}","page":"CellBase.jl","title":"CellBase.cellpar2mat","text":"convert cell parameters to column vectors\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.isvalidcellpar-NTuple{6, Any}","page":"CellBase.jl","title":"CellBase.isvalidcellpar","text":"Check if cell parameters are valid\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.vec2cellpar-Tuple{AbstractVector{T} where T, AbstractVector{T} where T, AbstractVector{T} where T}","page":"CellBase.jl","title":"CellBase.vec2cellpar","text":"Convert cell vectors to cell parameters\n\nReturns an static array of the cell parameters\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.volume-NTuple{6, Any}","page":"CellBase.jl","title":"CellBase.volume","text":"Compute volume from cell parameters\n\n\n\n\n\n","category":"method"},{"location":"#Neighbour-lists","page":"CellBase.jl","title":"Neighbour lists","text":"","category":"section"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"Type and functions for handling neighbour lists.","category":"page"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"Modules = [CellBase]\nOrder = [:type, :function]\nPages = [\"neighbour.jl\"]","category":"page"},{"location":"#CellBase.ExtendedPointArray","page":"CellBase.jl","title":"CellBase.ExtendedPointArray","text":"Represent an array of points after expansion by periodic boundary\n\n\n\n\n\n","category":"type"},{"location":"#CellBase.ExtendedPointArray-Tuple{Cell, Any}","page":"CellBase.jl","title":"CellBase.ExtendedPointArray","text":"ExtendedPointArray(cell::Cell, rcut)\n\nConstructed an ExtendedPointArray from a given structure\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.NLIterator","page":"CellBase.jl","title":"CellBase.NLIterator","text":"Iterator interface for going through all neighbours\n\n\n\n\n\n","category":"type"},{"location":"#CellBase.NLIteratorWithVector","page":"CellBase.jl","title":"CellBase.NLIteratorWithVector","text":"Iterator interface for going through all neighbours\n\n\n\n\n\n","category":"type"},{"location":"#CellBase.NeighbourList","page":"CellBase.jl","title":"CellBase.NeighbourList","text":"Type for representing a neighbour list\n\n\n\n\n\n","category":"type"},{"location":"#CellBase.NeighbourList-2","page":"CellBase.jl","title":"CellBase.NeighbourList","text":"NeighbourList(ea::ExtendedPointArray, rcut, nmax=100; savevec=false)\n\nConstruct a NeighbourList from an extended point array for the points in the original cell\n\n\n\n\n\n","category":"type"},{"location":"#CellBase.eachneighbour-Tuple{NeighbourList, Any}","page":"CellBase.jl","title":"CellBase.eachneighbour","text":"Iterate the neighbours of a site in the original cell. Returns a tuple of (originalindex, extendedindex, distance) for each iteration\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.eachneighbourvector-Tuple{NeighbourList, Any}","page":"CellBase.jl","title":"CellBase.eachneighbourvector","text":"Iterate the neighbours of a site in the original cell. Returns a tuple of (originalindex, extendedindex, distance, vector) for each iteration\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.nions_extended-Tuple{ExtendedPointArray}","page":"CellBase.jl","title":"CellBase.nions_extended","text":"Number of ions in the extended cell\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.nions_extended-Tuple{NeighbourList}","page":"CellBase.jl","title":"CellBase.nions_extended","text":"Number of ions in the extended cell\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.nions_orig-Tuple{ExtendedPointArray}","page":"CellBase.jl","title":"CellBase.nions_orig","text":"Number of ions in the original cell\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.nions_orig-Tuple{NeighbourList}","page":"CellBase.jl","title":"CellBase.nions_orig","text":"Number of ions in the original cell\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.num_neighbours-Tuple{NeighbourList, Any}","page":"CellBase.jl","title":"CellBase.num_neighbours","text":"Number of neighbours for a point\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.rebuild!-Tuple{ExtendedPointArray, Any}","page":"CellBase.jl","title":"CellBase.rebuild!","text":"Rebuild the ExtendedPointArray for an existing cell\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.rebuild!-Tuple{NeighbourList, Cell}","page":"CellBase.jl","title":"CellBase.rebuild!","text":"rebuild(nl::NeighbourList, cell::Cell)\n\nPerform a full rebuild of the NeighbourList with the latest geometry of the cell\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.rebuild!-Tuple{NeighbourList, ExtendedPointArray}","page":"CellBase.jl","title":"CellBase.rebuild!","text":"rebuild!(nl::NeighbourList, ea::ExtendedPointArray)\n\nPerform a full rebuild of the neighbour list from scratch\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.update!-Tuple{NeighbourList, Cell}","page":"CellBase.jl","title":"CellBase.update!","text":"update!(nl::NeighbourList, cell::Cell)\n\nUpdate the NeighbourList with the latest geometry of the Cell. No rebuilding is performed\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.update!-Tuple{NeighbourList}","page":"CellBase.jl","title":"CellBase.update!","text":"Update the calculated distances/vectors but do not rebuild the whole neighbour list\n\n\n\n\n\n","category":"method"},{"location":"#IO","page":"CellBase.jl","title":"IO","text":"","category":"section"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"Code for read/writing file of crystal structures.","category":"page"},{"location":"","page":"CellBase.jl","title":"CellBase.jl","text":"Modules = [CellBase, CellBase.CellIO, CellBase.DotCastep, CellBase.SheapIO]\nOrder = [:type, :function]\nPages = [\"io/io_cell.jl\", \"io/io_res.jl\", \"io/io_xyz.jl\", \"io/io_dotcastep.jl\", \"io/io_sheap.jl\"]","category":"page"},{"location":"#CellBase.read_res-Tuple{Vector{String}}","page":"CellBase.jl","title":"CellBase.read_res","text":"Read an array containing the lines of the SHELX file\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.write_res-Tuple{AbstractString, Cell}","page":"CellBase.jl","title":"CellBase.write_res","text":"write_res(fname::AbstractString, structure::Cell)\n\nWrite out SHELX format data to a file\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.write_res-Tuple{IO, Cell}","page":"CellBase.jl","title":"CellBase.write_res","text":"write_res(io::IO, structure::Cell)\n\nWrite out SHELX format data\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.push_xyz!-Tuple{Any, Cell}","page":"CellBase.jl","title":"CellBase.push_xyz!","text":"xyz lines for a single frame\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.write_xyz-Union{Tuple{T}, Tuple{Any, Array{Cell{T}, 1}}} where T","page":"CellBase.jl","title":"CellBase.write_xyz","text":"Write snapshots to a xyz file\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.CellIO._parse_tagline-Tuple{Any, Any}","page":"CellBase.jl","title":"CellBase.CellIO._parse_tagline","text":"Parse the tag line, return a dictionary of the tagline is in the form of <ionsetname> % KEY=VALUE KEY=VALUE\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.CellIO.clean_lines-Tuple{Any}","page":"CellBase.jl","title":"CellBase.CellIO.clean_lines","text":"Clean lines by removing new line symbols and comments\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.CellIO.filter_block-Tuple{Any}","page":"CellBase.jl","title":"CellBase.CellIO.filter_block","text":"Get rid of blocks \n\n\n\n\n\n","category":"method"},{"location":"#CellBase.CellIO.filter_block-Union{Tuple{T}, Tuple{Any, Vector{T}}} where T<:AbstractString","page":"CellBase.jl","title":"CellBase.CellIO.filter_block","text":"Filter away block names contained in only\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.CellIO.find_block-Tuple{Any, Any}","page":"CellBase.jl","title":"CellBase.CellIO.find_block","text":"Find the block with given name, return a Vector of the lines Comments are skipped\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.CellIO.parse_taglines-Union{Tuple{Vector{T}}, Tuple{T}} where T<:AbstractString","page":"CellBase.jl","title":"CellBase.CellIO.parse_taglines","text":"Parse tag lines, return a vector of dictionary containing parsed tags \n\n\n\n\n\n","category":"method"},{"location":"#CellBase.CellIO.read_cell-Tuple{AbstractString}","page":"CellBase.jl","title":"CellBase.CellIO.read_cell","text":"Read a cell file\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.CellIO.read_cell-Union{Tuple{Vector{T}}, Tuple{T}} where T<:AbstractString","page":"CellBase.jl","title":"CellBase.CellIO.read_cell","text":"Read content of a cell file\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.CellIO.read_cellmat-Tuple{Any}","page":"CellBase.jl","title":"CellBase.CellIO.read_cellmat","text":"Read cell related sections\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.CellIO.read_num_block-Union{Tuple{T}, Tuple{Vector{T}, Int64}} where T<:AbstractString","page":"CellBase.jl","title":"CellBase.CellIO.read_num_block","text":"Read a numerical block in the form of a vector of String\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.CellIO.read_positions-Tuple{Any, Any}","page":"CellBase.jl","title":"CellBase.CellIO.read_positions","text":"Read positions related sections\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.CellIO.read_seed-Tuple{AbstractString}","page":"CellBase.jl","title":"CellBase.CellIO.read_seed","text":"Read in cell \n\n\n\n\n\n","category":"method"},{"location":"#CellBase.CellIO.read_seed-Union{Tuple{Vector{T}}, Tuple{T}} where T<:AbstractString","page":"CellBase.jl","title":"CellBase.CellIO.read_seed","text":"Read seed with the label parsed. At this stage no explansions are done\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.DotCastep.SnapShot","page":"CellBase.jl","title":"CellBase.DotCastep.SnapShot","text":"Representation of a snapshot\n\n\n\n\n\n","category":"type"},{"location":"#CellBase.DotCastep._get_nions-Tuple{Any}","page":"CellBase.jl","title":"CellBase.DotCastep._get_nions","text":"Get the number of ions\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.DotCastep._read_snapshot","page":"CellBase.jl","title":"CellBase.DotCastep._read_snapshot","text":"Read the initial structure and its energy/forces\n\n\n\n\n\n","category":"function"},{"location":"#CellBase.DotCastep.read_castep-Tuple{String}","page":"CellBase.jl","title":"CellBase.DotCastep.read_castep","text":"read_castep(fname::String)\n\nRead a CASTEP file, return a Vector of the Snapshots\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.DotCastep.read_coord_table!","page":"CellBase.jl","title":"CellBase.DotCastep.read_coord_table!","text":"Read fractional coordinates (column vectors) from a table\n\n\n\n\n\n","category":"function"},{"location":"#CellBase.DotCastep.read_lattice-Tuple{Any}","page":"CellBase.jl","title":"CellBase.DotCastep.read_lattice","text":"Read lattice vectors from the lines (column vectors)\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.DotCastep.skip_to_header-Tuple{Vector{String}, Int64, Int64}","page":"CellBase.jl","title":"CellBase.DotCastep.skip_to_header","text":"Skip to the next header\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.SheapIO.SheapOptions","page":"CellBase.jl","title":"CellBase.SheapIO.SheapOptions","text":"Commandline options for SHEAP\n\n\n\n\n\n","category":"type"},{"location":"#CellBase.SheapIO._SheapDefaultOptions","page":"CellBase.jl","title":"CellBase.SheapIO._SheapDefaultOptions","text":"Record the default options of SHEAP commandline\n\n\n\n\n\n","category":"type"},{"location":"#CellBase.SheapIO.parse_sheap_output-Tuple{IO}","page":"CellBase.jl","title":"CellBase.SheapIO.parse_sheap_output","text":"Prase the output of SHEAP from an IO object\n\n\n\n\n\n","category":"method"},{"location":"#CellBase.SheapIO.run_sheap-Tuple{Any, CellBase.SheapIO.SheapOptions}","page":"CellBase.jl","title":"CellBase.SheapIO.run_sheap","text":"run_sheap(vecs, opt::SheapOptions;metadata=repeat([SheapMetadata()], length(vecs)), show_stderr=true)\n\nRun SHEAP for a iterator of vectors with the given options.\n\n\n\n\n\n","category":"method"}]
}
