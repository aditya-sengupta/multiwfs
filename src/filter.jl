using Base: prod
using StaticArrays

abstract type Filter end

function f2s(f)
    return 1im * 2.0 * π * f
end

function ar1_coeff(f_cutoff, f_loop)
    return exp(-2π * f_cutoff / f_loop)
end

function ar1(f_cutoff, f_loop, filter_type)
    α = ar1_coeff(f_cutoff, f_loop)
    if filter_type == "high"
        return ZPKFilter(1, α, α)
    elseif filter_type == "low"
        return ZPKFilter(0, α, 1 - α)
    end
end

mutable struct ZPKFilter <: Filter
    z::Float64
    p::ComplexF64
    k::Float64
    x_nm1::ComplexF64
    y_nm1::ComplexF64

    function ZPKFilter(z, p, k)
        new(z, p, k, 0.0, 0.0)
    end
end

function output!(zpkf::ZPKFilter, x_n)
    y_n = zpkf.p * zpkf.y_nm1 + zpkf.k * x_n - zpkf.k * zpkf.z * zpkf.x_nm1
    zpkf.x_nm1 = x_n
    zpkf.y_nm1 = y_n
    return y_n
end

function reset!(zpkf::ZPKFilter)
    zpkf.x_nm1 = 0.0
    zpkf.y_nm1 = 0.0
end

function transfer_function(zpkf::ZPKFilter, f_over_f_loop)
    z = exp(2π * im * f_over_f_loop)
    return zpkf.k * (z - zpkf.z) / (z - zpkf.p)
end

struct CascadeZPKFilter <: Filter
    filters::Vector{ZPKFilter}

    function CascadeZPKFilter(z::Vector{Float64}, p::Vector{ComplexF64}, k::Float64)
        filters = []
        k_i = 1
        for (i, (z_i, p_i)) in enumerate(zip(z, p))
            if i == length(z)
                k_i = k
            end
            push!(filters, ZPKFilter(z_i, p_i, k_i))
        end
        new(filters)
    end
end

function output!(czpkf::CascadeZPKFilter, x)
    for f in czpkf.filters
        x = output!(f, x)
    end
    return x
end

function transfer_function(czpkf::CascadeZPKFilter, f_over_f_loop)
    return mapreduce(filt -> transfer_function(filt, f_over_f_loop), *, czpkf.filters)
end

function reset!(czpkf::CascadeZPKFilter)
    for f in czpkf.filters
        reset!(f)
    end
end

export ar1, ZPKFilter, CascadeZPKFilter, transfer_function, output!, reset!
