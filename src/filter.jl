using Plots

function α(f_cutoff, f_loop)
    return exp(-2π * f_cutoff / f_loop)
end

mutable struct ZPKFilter
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

struct CascadeZPKFilter
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

function reset!(czpkf::CascadeZPKFilter)
    for f in czpkf.filters
        f.x_nm1 = 0.0
        f.y_nm1 = 0.0
    end
end

czpkf = CascadeZPKFilter([1.0, 1.0, 1.0, 1.0], ComplexF64[0.9803579415531668 - 0.08785432565973045im, 0.8896705374015489 - 0.09630554662137608im, 0.8896705374015489 + 0.09630554662137608im, 0.9803579415531668 + 0.08785432565973045im], 0.8353022013649711)

begin
    outs = []
    reset!(czpkf)
    freq = 3.0
    times = 0.0:1e-2:(2/freq * 10)
    ins = sin.(2π .* freq .* times)
    for x in ins
        push!(outs, output!(czpkf, x))
    end

    plot(ins)
    plot!(real.(outs))
end