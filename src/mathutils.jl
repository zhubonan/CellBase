#=
Utility modules - contains common routines
=#
using LinearAlgebra
   
const dgrd = π / 180.
const rddg = 1 / dgrd

function angle(x::AbstractVector, y::AbstractVector, radian=false)
    lx = norm(x)
    ly = norm(y)
    dotp = dot(x, y)
    cos_ang = dotp / lx / ly
    ang = acos(cos_ang) 
    if !radian
        ang *= 180.0 / pi
    end
    return ang
end


"""
Convert cell vectors to cell parameters

Returns an static array of the cell parameters
"""
function vec2cellpar(va::AbstractVector, vb::AbstractVector, vc::AbstractVector)
    a = norm(va)
    b = norm(vb)
    c = norm(vc)

    α = angle(vb, vc)
    β = angle(va, vc)
    γ = angle(va, vb)
    return [a, b, c, α, β, γ]
end

function vec2cellpar(cell::AbstractMatrix)
    va = cell[:, 1]
    vb = cell[:, 2]
    vc = cell[:, 3]
    return vec2cellpar(va, vb, vc)
end

 
function _cellpar_trig(α, β, γ)
    # For the orthorhombic cell case
    eps = 1.4e-14   # Error tolorance
    # Alpha
    if abs(abs(α) - 0.5π) < eps
        cos_alpha = 0.0
    else
        cos_alpha = cos(α)
    end
    # β
    if abs(abs(β) - 0.5π) < eps
        cos_beta = 0.0
    else
        cos_beta = cos(β)
    end

    if abs(γ - 0.5π) < eps
        cos_gamma = 0.0
        sin_gamma = 1.0
    elseif abs(γ + 0.5π) < eps
        cos_gamma = 0.0
        sin_gamma = -1.0
    else
        cos_gamma = cos(γ)
        sin_gamma = sin(γ)
    end
    return cos_alpha, cos_beta, cos_gamma, sin_gamma
end

radian(a, b, c) = dgrd * a, dgrd * b, dgrd * c  

"Check if cell parameters are valid"
function isvalidcellpar(a, b, c, α, β, γ; degree=true)
    if degree
        α, β, γ = radian(α, β, γ)
    end
    cos_alpha, cos_beta, cos_gamma, sin_gamma = _cellpar_trig(α, β, γ)
    cx = cos_beta
    cy = (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    tmp = 1 - cx * cx - cy * cy
    tmp >= 0.
end


"convert cell parameters to column vectors"
function cellpar2mat(a, b, c, α, β, γ; degree=true)

    if degree
        α, β, γ = radian(α, β, γ)
    end
    a_direction = [1., 0., 0.]
    ab_normal = [0., 0., 1.]

    cos_alpha, cos_beta, cos_gamma, sin_gamma = _cellpar_trig(α, β, γ)

    # Buildce cell vectors
    va = [a, 0., 0.]
    vb = [cos_gamma * b, sin_gamma * b, 0.0]
    cx = cos_beta
    cy = (cos_alpha - cos_beta * cos_gamma) / sin_gamma 
    tmp = 1. - cx * cx - cy * cy
    @assert tmp > 0 "Cell parameters are not valid"
    cz = sqrt(tmp)
    vc = [cx * c, cy * c, cz * c]

    return hcat(va, vb, vc)
end

"Compute volume from cell parameters"
function volume(a, b, c, α, β, γ; degree=true)
    if degree
        α, β, γ = radian(α, β, γ)
    end
    a * b * c * sqrt(1 + 2 * cos(α) * cos(β) * cos(γ) - cos(α)^2 - cos(β)^2 - cos(γ)^2)
end
