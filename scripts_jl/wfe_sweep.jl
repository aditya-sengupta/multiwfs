using multiwfs
using multiwfs: Hrej
using Plots
using ProgressMeter
using Base.Threads
using QuadGK

sys = AOSystem(1000.0, 1.0, 0.01, 0.999, 10, "high", 1.0)

gains = gain_map(sys)

function wfe(sys, fmin=0.1, fmax=500.0)
    return quadgk(
        f -> abs2(Hrej(sys, f)) * f^(-11/3),
        fmin, fmax
    )[1]
end

function wfe_sweep(sys, gains; f_cutoffs = 0.1:0.1:100.0, delays = 0.0:0.1:1.0)
    wfes = zeros(length(f_cutoffs), length(delays));
    @showprogress @threads for (i, fc) in collect(enumerate(f_cutoffs))
        for (j, d) in enumerate(delays)
            tsys = AOSystem(sys.f_loop, d, gains[i,j], sys.leak, sys.fpf, sys.filter_type, fc)
            wfes[i,j] = wfe(tsys)
        end
    end
    f_cutoffs, delays, wfes
end

function plot_wfe_sweep(f_cutoffs, delays, wfes; t="", kwargs...)
    p = plot(xlabel="Cutoff frequency (Hz)", ylabel="Relative WFE", title=t; kwargs...)
    c = palette([:blue, :red], length(delays))
    for (i, (r, d)) in enumerate(zip(eachcol(wfes), delays))
        plot!(f_cutoffs, r, label="$d frames", c=c[i])
    end
    p
end

sys.filter_type = "high"
f_cutoffs, delays, wfes_high = wfe_sweep(sys, gains)
sys.filter_type = "low"
f_cutoffs, delays, wfes_low = wfe_sweep(sys, gains)

plot(
    plot_wfe_sweep(f_cutoffs, delays, wfes_low, t="LPF"; yscale=:log10),
    plot_wfe_sweep(f_cutoffs, delays, wfes_high, t="HPF"; legend=:bottomright),
)

plot_wfe_sweep(f_cutoffs, delays, max.(wfes_low, wfes_high), legend=:bottomright)