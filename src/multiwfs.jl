module multiwfs
    using ControlSystems
    using Plots

    function f2s(f)
        return 1im * 2.0 * π * f
    end

    mutable struct AOSystem
        f_loop::Float64
        Ts::Float64
        frame_delay::Float64
        τ::Float64
        gain::Float64
        leak::Float64
        fpf::Int64
        filter_type::String
        filter_cutoff::Float64
        α::Float64

        function AOSystem(f_loop, frame_delay, gain, leak, fpf, filter_type, filter_cutoff)
            Ts = 1 / f_loop
            τ = frame_delay * Ts
            frame_delay =  floor(frame_delay) + round((frame_delay-floor(frame_delay))*fpf)/fpf
            α = exp(-2π * filter_cutoff / f_loop)
            new(f_loop, Ts, frame_delay, τ, gain, leak, fpf, filter_type, filter_cutoff, α)
        end
    end

    function Hwfs(ao::AOSystem, s)
        return (1.0 - exp(-ao.Ts * s)) / (ao.Ts * s)
    end

    function Hzoh(ao::AOSystem, s)
        return Hwfs(ao, s)
    end

    function Hlag(ao::AOSystem, s)
        return exp(-ao.τ * s)
    end

    function Hint(ao::AOSystem, s)
        return ao.gain / (1.0 - exp(-ao.Ts * s)) # integrator
    end

    function Hlint(ao::AOSystem, s)
        return ao.gain / (1.0 - ao.leak * exp(-ao.Ts * s)) # leaky integrator
    end

    function Hcont(ao::AOSystem, s)
        return Hlint(ao, s)
    end

    function Hfilter(ao::AOSystem, s)
        if ao.filter_type == "high"
            return ao.α * (1 - exp(-s / ao.f_loop)) / (1 - ao.α * exp(-s / ao.f_loop))
        elseif ao.filter_type == "low"
            return (1 - ao.α) / (1 - ao.α * exp(-s / ao.f_loop))
        else
            return 1
        end
    end

    function Hol(ao::AOSystem, s::TransferFunction)
        return Hwfs(ao, s) * Hlag(ao, s) * Hcont(ao, s) * Hfilter(ao, s) * Hzoh(ao, s)
    end

    function Hol(ao::AOSystem, f)
        s = f2s(f)
        return Hwfs(ao, s) * Hlag(ao, s) * Hcont(ao, s) * Hfilter(ao, s) * Hzoh(ao, s)
    end

    function Hrej(ao::AOSystem, f)
        return 1.0 / (1.0 + Hol(ao, f))
    end

    function Hcl(ao::AOSystem, f)
        Hol_val = Hol(ao, f)
        return Hol_val / (1.0 + Hol_val)
    end

    function Hn(ao::AOSystem, f)
        return Hcl(ao, f) / Hwfs(ao, f2s(f))
    end

    unfiltered = AOSystem(1000.0, 0.5, 0.2, 0.99, 10, "none", 0.0)
    low_sys = AOSystem(1000.0, 0.5, 0.2, 0.99, 10, "low", 10.0)
    high_sys = AOSystem(1000.0, 0.5, 0.2, 0.99, 10, "high", 10.0)

end # module multiwfs
