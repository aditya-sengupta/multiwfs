using multiwfs
using Plots
using DSP
import multiwfs: Hrej, Hol, Hcont, Hfilter, Hwfs, Hlag, Hzoh, Hlint

function Hrej(systems::Vector{AOSystem{T}}, f) where T
    ol = sum(Hol(sys, f) for sys in systems)
    return 1 / (1 + ol)
end

f_loop = 1000.0
f_cutoff = 50.0
ar1_low = ar1_filter(f_cutoff, f_loop, "low")
ar1_high = ar1_filter(f_cutoff, f_loop, "high")
# parameters are:
# f_loop, frame_delay, gain, leak, fpf
sys_low = AOSystem(f_loop, 1.0, 1.3, 0.999, 1, ar1_low)
sys_high = AOSystem(f_loop, 1.0, 0.61, 0.999, 1, ar1_high)

f = 0.032:0.032:500.0
p1 = begin
    plot(legend=:bottomright, xticks=[1e-1, 1e0, 1e1, 1e2], yticks=[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0], ylabel="|Rejection transfer function|²")
    plot_psd_p!(f, abs2.(Hrej.(Ref(sys_low), f)), label="LPF", color=:purple)
    plot_psd_p!(f, abs2.(Hrej.(Ref(sys_high), f)), label="HPF", color=:royalblue4)
    plot_psd_p!(f, abs2.(Hrej.(Ref([sys_high, sys_low]), f)), label="LPF+HPF", color=:cadetblue)
end

# what I want here is for phi_to_X to match up with Hrej on both filters
# if I plug in the appropriate controller TFs
# I'll play around with z^-1 terms to make this match

begin
    s = 2π * im * f
    phi_to_X.(s / f_loop,
        sT -> Hlint(sys_high, sT * f_loop) * Hfilter(sys_high, sT * f_loop),
        sT -> Hlint(sys_low, sT * f_loop) * Hfilter(sys_low, sT * f_loop),
        1
    ) |> p -> plot(f, abs2.(p), label="two-WFS plant", color=:red, legend=:bottomleft)
    plot_psd_p!(f, abs2.(1 ./ (1 .+ (Hol.(Ref(sys_high), f) .+ Hol.(Ref(sys_low), f)))), label="LPF+HPF", color=:cadetblue, ylabel="ETF")
end

