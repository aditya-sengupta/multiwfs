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
    Ts::Float64

    function LQG(A, D, C, K, L, Ts)
        ikcA = (I - K * C) * A
        ikcD = (I - K * C) * D
        new(A, D, C, K, L, ikcA, ikcD, Ts)
    end
end

function transfer_function(lqg::LQG, s::Complex)
    zinv = exp(-s*lqg.Ts)
    internal_term = inv(I - lqg.ikcA * zinv)
    numerator = (lqg.L * internal_term * lqg.K)[1,1]
	denominator = (I - lqg.L * internal_term * lqg.ikcD * zinv)[1,1]
    return numerator / denominator
end

struct FilteredLQG{F} <: Controller
    lqg::LQG
    cfilter::F
end

function transfer_function(fl::FilteredLQG, s::Complex)
    return transfer_function(fl.lqg, s) * transfer_function(fl.cfilter, fl.lqg.Ts * s)
end

export FilteredIntegrator, LQG, FilteredLQG, transfer_function