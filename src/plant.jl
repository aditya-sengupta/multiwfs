mutable struct Plant
    f_loop::Float64
    Ts::Float64
    frame_delay::Float64
    τ::Float64
    fpf::Int64

    function Plant(f_loop, frame_delay, fpf)
        Ts = 1 / f_loop
        frame_delay = floor(frame_delay) + round((frame_delay-floor(frame_delay))*fpf)/fpf
        τ = frame_delay * Ts
        new(f_loop, Ts, frame_delay, τ, fpf)
    end
end

function Hwfs(p::Plant, s)
    return (1.0 - exp(-p.Ts * s)) / (p.Ts * s)
end

function Hzoh(p::Plant, s)
    return Hwfs(p, s)
end

function Hlag(p::Plant, s)
    return exp(-p.τ * s)
end

function transfer_function(p::Plant, s::Complex)
    return Hwfs(p, s) * Hlag(p, s) * Hzoh(p, s)
end

export Plant