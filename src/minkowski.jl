#=
Perform Minkowski reduction of the given lattice basis (3 dimensional)
=#
using LinearAlgebra


function minkowski_reduce3d(B::AbstractMatrix)
    H = Float64[
       1 0 0
       0 1 0
       0 0 1
       ] 
    norms = map(norm, eachcol(B))
    MAX_IT = 10000
    for it in 1:MAX_IT
        # Sort vectors by norm
        H = H[:, sortperm(norms)]

        # Gauess-reduce smallest two vectors
        hw = H[:, 3]
        hu, hv = reduction_gauss(B, H[:, 1], H[:, 2])
        H = hcat(hu, hv, hw)
        R = B * H

        # Orthogonalize vectors using Gram-Schmidt
        u = R[:, 1]
        v = R[:, 2]
        X = u ./ norm(u)
        Y = v .- X * dot(v, X)
        Y ./= norm(Y)

        # Find cloest vector to last element of R
        tmp = vcat(rowvec(X), rowvec(Y)) * R
        pu, pv, pw = tmp[:, 1], tmp[:, 2], tmp[:, 3]
        nb = closest_vector(pw, pu, pv)

        # Update the basis
        H[:, 3] = H * [nb[1], nb[2], 1]
        R = B * H

        norms = map(norm, eachcol(R))
        if norms[3] >= norms[2]
            return R, H
        end 
    end
    throw(ErrorException("Cannot find reduced basis after  $(MAX_IT) iterations"))
end


"""
2D Gauss reduction

Works for column vectors instead of row vectors
"""
function reduction_gauss(B::Matrix, hu::Vector, hv::Vector)
    u = B * hu
    v = B * hv
    MAX_IT = 10000
    for it in 1:MAX_IT
        x = Int(round(dot(u, v) / dot(u, u)))
        hu_bak = copy(hu)
        hu .= hv .- x * hu
        hv = hu_bak
        u = B * hu
        v = B * hv
        if dot(u, u) >= dot(v, v)
            return hv, hu
        end
    end
    throw(ErrorException("Cannot find reduced basis after  $(MAX_IT) iterations"))
end
rowvec(vec::Vector) = collect(transpose(vec))

function closest_vector(t0, u, v)
    t = copy(t0)
    a = Int[0, 0]
    rs, cs = relavant_vectors_2d(u, v)
    dprev = Inf
    MAX_IT = 100000
    for it in 1:MAX_IT
        ds = map(norm, eachslice(rs .+ t, dims=2))
        index = argmin(ds)
        if index == 1 || ds[index] >= dprev
            return a
        end

        dprev = ds[index]
        r = rs[:, index]
        kopt = Int(round(-dot(t, r)/dot(r,r)))
        a .+= kopt * cs[:, index]
        t .= t0 .+ a[1] .* u + a[2] .* v
    end
    throw(ErrorException("Cannot find reduced basis after  $(MAX_IT) iterations"))
end


function relavant_vectors_2d(u, v)
    cs = [
        -1 -1 -1  0  0  0  1  1  1
        -1  0  1 -1  0  1 -1  0  1
    ]
    vs = hcat(u, v) * cs
    indices = sortperm(map(norm, eachcol(vs)))[1:7]
    return vs[:, indices], cs[:, indices]
end