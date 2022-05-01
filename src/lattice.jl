using Printf

export Lattice, mic
"""
Lattice{T}

The type represents a lattice in 3D space formed by three column lattice vectors.
"""
struct Lattice{T}
    matrix::Matrix{T}
    rec::Matrix{T}
end

function Lattice(matrix::Matrix{T}) where T
    Lattice(matrix, inv(matrix))
end

function Lattice(matrix::Matrix{T}) where {T<:Integer}
    Lattice(convert.(Float64, matrix), inv(matrix))
end

"""Construct from lattice vectors"""
function Lattice(va::Vector{T}, vb::Vector{T}, vc::Vector{T}) where T
    cellmat = hcat(va, vb, vc)
    Lattice(cellmat)
end

"""Construct from lattice parameters"""
function Lattice(cellpar::Vector{T}) where T
    cellmat = cellpar2mat(cellpar...)
    Lattice(cellmat)
end

"""Construct from lattice parameters"""
function Lattice(a::T, b::T, c::T, α::T, β::T, γ::T) where T
    cellmat = cellpar2mat(a, b, c, α, β, γ)
    Lattice(cellmat)
end

"""Construct from three lattice parameters"""
function Lattice(a::T, b::T, c::T) where T <: Real
    cellmat = zeros(T, 3, 3)
    cellmat[1, 1] = a
    cellmat[2, 2] = b
    cellmat[3, 3] = c
    Lattice(cellmat)
end

"Reciprocal lattice (no 2pi factor)"
reciprocal(l::Lattice) = l.rec

"Reciprocal lattice (no 2pi factor)"
rec_cellmat(l::Lattice) = l.rec

"""Get a matrix of column vectors"""
cellmat(lattice::Lattice) = lattice.matrix

"""Get a matrix of column vectors"""
cellmat_row(lattice::Lattice) = copy(transpose(lattice.matrix))

"""Lattice vectors"""
cellvecs(lattice::Lattice) = lattice.matrix[:, 1], lattice.matrix[:, 2], lattice.matrix[:, 3]

"""Volume of the cell"""
function volume(lattice::Lattice)
    a, b, c = cellvecs(lattice)
    dot(a, cross(b, c))
end

"""Cell parameters"""
cellpar(lattice::Lattice) = vec2cellpar(cellvecs(lattice)...)

"Fraction positions of a site"
frac_pos(s::Site, l::Lattice) = l.rec * s.position 

"Wrap a site back into the box"
function wrap!(s::Site, l::Lattice)
    frac = frac_pos(s, l)
    for i in 1:3
        v = frac[i] % 1
        v < 0 && (v += 1)
        frac[i] = v
    end 
    s.position[:] = l.matrix * frac
end

function Base.show(io::IO, l::Lattice)
    println("Lattice: ")
    cellmat = l.matrix
    for (i, s) in zip(1:3, ['a', 'b', 'c'])
        @printf "%s %5.3f   %5.3f   %5.3f\n" s  cellmat[1, i] cellmat[2, i] cellmat[3, i]
    end
end


"""
    random_vec_in_cell(cell::Matrix{T}; scales::Vector=ones(size(cell)[1])) where T

Get an randmo vector within a unit cell. The cell is a matrix made of column vectors.
"""
function random_vec_in_cell(cell::Matrix{T}; scales::Vector=ones(size(cell)[1])) where T
    out = zeros(T,  size(cell)[2])
    for (i, scale) in enumerate(scales)
        out[:] += cell[:, i] * (scale * randf())
    end
    out
end

"Random displacement vector in a unit cell"
function random_vec(l::Lattice)
    random_vec_in_cell(cellmat(l))
end

"Return scaled positions"
get_scaled_positions(l::Lattice, v::AbstractVecOrMat) = reciprocal(l) * v

"Alias  to `get_scaled_positions`"
const scaled_positions = get_scaled_positions

function wrap_positions(l::Lattice, v::AbstractVecOrMat)
    scaled = get_scaled_positions(l, v)
    scaled .-= floor.(scaled)
    cellmat(l) * scaled
end


"""
Compute naive MIC representation of the vector(s)

Not safe for skewed cells. Requires

```julia
norm(length) < 0.5 * min(cellpar(l)[1:3])
```

See also:
- W. Smith, "The Minimum Image Convention in Non-Cubic MD Cells", 1989,
    http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.57.1696.

"""
function mic_naive(l::Lattice, v::AbstractMatrix)
    scaled = scaled_positions(l, v)
    scaled .-= floor.(scaled .+ 0.5)
    vmin = cellmat(l) * scaled
    vlen = map(norm, eachcol(vmin))
    vmin, vlen
end


"""
     mic_safe(l::Lattice, v::AbstractVector)

Compute MIC representation of vectors using a safe approach based on Minkowski reduced cell
"""
mic_naive(l::Lattice, v::AbstractVector) = mic_naive(l, reshape(v, 3, 1))


"""
     mic_safe(l::Lattice, v::AbstractMatrix)

Compute MIC representation of vectors using a safe approach based on Minkowski reduced cell
"""
function mic_safe(l::Lattice, v::AbstractMatrix)
    rcell, _ = minkowski_reduce3d(cellmat(l))
    rl = Lattice(rcell)
    wrapped = wrap_positions(rl, v)


    shiftvec = mic_shiftvecs(rcell)
    minvecs = similar(wrapped)
    mic_dist = zeros(size(wrapped)[2])

    # Search for minimum distances
    for (i, vec) in enumerate(eachcol(wrapped))
        mind2 = dot(vec, vec)
        minshift = 1  # Default to one as we know the first shift is 0, 0, 0
        # Find the shorts vector
        for (ishift, svec) in enumerate(eachcol(shiftvec))
            d2 = 0.
            # dot(svec .+ vec, svec .+ vec)
            for ii = 1:3
                @inbounds tmp = svec[ii] + vec[ii]
                d2 += tmp * tmp
            end
            # Shorter equivalent vector detected - store the length and indices
            if  d2 < mind2
                # Store the shifted vecor
                minshift = ishift 
                mind2 = d2
            end
        end
        # Store the shifted vector and its length
        for ii = 1:3
            minvecs[ii, i] = vec[ii] + shiftvec[ii, minshift]
        end
        mic_dist[i] = sqrt(mind2)
    end
    minvecs, mic_dist
end

"""
     mic_safe(l::Lattice, v::AbstractVector)

Compute MIC representation of vectors using a safe approach based on Minkowski reduced cell
"""
mic_safe(l::Lattice, v::AbstractVector) = mic_safe(l, reshape(v, 3, 1))

function mic_shiftvecs(rcell::AbstractMatrix)
    # Assign shift vectors
    # For the reduced cell we just shift by (-1, 1) in each dimension
    shifts = zeros(Int, 3, 28)
    i=2
    for tup in Base.product([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])
        shifts[:, i] .= tup
        i +=1
    end
    rcell * shifts
end


function mic(l::Lattice, vec::AbstractVecOrMat)

    vmic, dmic = mic_naive(l, vec)
    a, b, c, _, _, _ = cellpar(l)
    thresh = min(a, b, c) / 2
    if any(x -> x > thresh, dmic)
        vmic, dmic = mic_safe(l, vec)
    end
    return vmic, dmic
end