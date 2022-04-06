#= 
functions for handling periodic conditions
=#
using LinearAlgebra

"""
    shift_vectors(lattice::AbstractMatrix, shift1, shift2, shift3)

Compute the shift vectors needed to include all parts of the lattice within a
cut off radius.
"""
function shift_vectors(lattice::AbstractMatrix, s1max, s2max, s3max, s1min=-s1max, s2min=-s2max, s3min=-s3max)::Matrix{Float64}
    a1, a2, a3, = lattice[:, 1], lattice[:, 2], lattice[:, 3]
    nshifts = (s1max - s1min + 1) * (s2max - s2min + 1) * (s3max - s3min +1)
    shift_vectors = Matrix{Float64}(undef, 3, nshifts)

    itmp = 1
    for s3 in s3min:s3max
        for s2 in s2min:s2max
            for s1 in s1min:s1max
                for i in 1:3
                    shift_vectors[i, itmp] = s1 * a1[i] + s2 * a2[i] + s3 * a3[i]
                end
                itmp += 1
            end  # s1
        end  # s2
    end  # s3
    return shift_vectors
end

"Compute shift vectors given lattice matrix and cut off radius"
function shift_vectors(lattice::AbstractMatrix, rc::Real; safe=false)::Matrix{Float64}
    shifts = max_shifts(lattice, rc; safe=safe)
    shift_vectors(lattice, shifts...)
end

"""
    shift_indices(lattice::AbstractArray, rc::Real)

Compute the shift indixing vectors:
  [1 1 1 ...
   1 1 1 ...
   1 2 3 ... 
  ]
"""
function shift_indices(lattice::AbstractArray, rc::Real)
    shift1max, shift2max, shift3max = max_shifts(lattice, rc)
    shift_indices(shift1max, shift2max, shift3max)
end

function shift_indices(shift1::Int, shift2::Int, shift3::Int)
    nshifts = (shift1 * 2 + 1) * (shift2 * 2 + 1) * (shift3 * 2 +1)
    idx = Array{Int}(undef, (3, nshifts))
    itmp = 1
    for s1 in -shift1:shift1, s2 in -shift2:shift2, s3 in -shift3:shift3
        idx[1, itmp] = s1
        idx[2, itmp] = s2
        idx[3, itmp] = s3
        itmp += 1
    end
    idx
end


"""
Compute the require periodicities required to fill a cut off sphere of a certain radius
"""
function max_shifts(lattice::AbstractArray, rc::Real; safe=false)

    if safe
        diag_vec = sum(lattice, dims=2)
        diag = sqrt(dot(diag_vec, diag_vec))
        rc += diag / 2
    end

    invl = inv(lattice)
    b1, b2, b3 = invl[:, 1], invl[:, 2], invl[:, 3]

    shift1max = ceil(Int, rc * norm(b1))
    shift2max = ceil(Int, rc * norm(b2))
    shift3max = ceil(Int, rc * norm(b3))
    return shift1max, shift2max, shift3max
end