using multiwfs
using NPZ
using CairoMakie

fastncp_ts = npzread("data/olt/olt_fastncp.npy")
k = 25102
std(fastncp[k:k+2048])
psd_l = gen_avg_per_unb(fastncp[k:k+2048], 2048)
sqrt(sum(psd_l))

function loglog(x, y=nothing; kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1], xscale=log10, yscale=log10; kwargs...)
    if isnothing(y)
        lines!(ax, x)
    else
        lines!(ax, x, y)
    end
    fig
end

loglog(range(0.0, 500.0, 1025)[2:end], psd_l[1:1024], xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)")