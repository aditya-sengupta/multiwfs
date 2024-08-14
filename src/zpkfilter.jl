using Base: prod
using StaticArrays

# much of this is already in DSP.jl
# but I'm re-implementing it for better control/understanding

function f2s(f)
    return 1im * 2.0 * π * f
end

struct ZPKFilter{N}
    z::SVector{N,Float64}
    p::SVector{N,ComplexF64}
    k::Float64
    prev_x::MVector{N,ComplexF64}
    prev_y::MVector{N,ComplexF64}

    function ZPKFilter(z::AbstractArray, p::AbstractArray, k::Number)
        N = length(z)
        x = @MVector zeros(ComplexF64,N)
        y = @MVector zeros(ComplexF64,N)
        new{N}(z, p, k, x, y)
    end

    function ZPKFilter(z::Number, p::Number, k::Number)
        ZPKFilter([z], [p], k)
    end
end

function output!(zpkf::ZPKFilter{N}, x_n) where N
    y_n = 0.0 + 0.0im # intermediate results can be complex but final ones shouldn't be!
    
    for i in 1:N
        k = (i == N ? zpkf.k : 1)
        y_n = zpkf.p[i] * zpkf.prev_y[i] + k * x_n - k * zpkf.z[i] * zpkf.prev_x[i]
        zpkf.prev_x[i] = x_n
        zpkf.prev_y[i] = y_n
        x_n = y_n
    end
    return real(y_n)
end

function reset!(zpkf::ZPKFilter{N}) where N
    zpkf.prev_x[:] = zeros(N)
    zpkf.prev_y[:] = zeros(N)
end

function transfer_function(zpkf::ZPKFilter, s)
    z = exp(s)
    return zpkf.k * prod((z - zv) / (z - p) for (zv, p) in zip(zpkf.z, zpkf.p))
end

function ar1_coeff(f_cutoff, f_loop)
    return exp(-2π * f_cutoff / f_loop)
end

function ar1_filter(f_cutoff, f_loop, filter_type)
    α = ar1_coeff(f_cutoff, f_loop)
    if filter_type == "high"
        return ZPKFilter(1, α, α)
    elseif filter_type == "low"
        return ZPKFilter(0, α, 1 - α)
    end
end

export ar1_filter, ZPKFilter, transfer_function, output!, reset!, f2s
