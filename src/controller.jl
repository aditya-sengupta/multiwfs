using LinearAlgebra: I

abstract type Controller end

include("make_lqg.jl")

mutable struct FilteredIntegrator{D,F} <: Controller
    gain::D
    leak::D
    cfilter::F
    Ts::Float64
end

function transfer_function(fi::FilteredIntegrator, s)
    fi.gain / (1.0 - fi.leak * exp(-fi.Ts * s)) * transfer_function(fi.cfilter, fi.Ts * s)
end

struct LQG <: Controller
    A::Matrix{Float64}
    D::Matrix{Float64}
    C::Matrix{Float64}
    K::Matrix{Float64}
    L::Matrix{Float64}
    ikcA::Matrix{Float64}
    ikcD::Matrix{Float64}

    function LQG(A, D, C, K, L)
        ikcA = (I - K * C) * A
        ikcD = (I - K * C) * D
        new(A, D, C, K, L, ikcA, ikcD)
    end
end

function transfer_function(lqg::LQG, s::Complex)
    zinv = exp(-s)
    numerator = (lqg.L * inv(I - lqg.ikcA * zinv))[1,1]
	denominator = (I - lqg.L * inv(I - lqg.ikcA * zinv) * lqg.ikcD * zinv)[1,1]
    return numerator * (zinv) / denominator
end

export FilteredIntegrator, LQG, transfer_function