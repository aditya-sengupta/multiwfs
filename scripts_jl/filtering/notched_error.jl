using multiwfs
using multiwfs: Hcont, Hfilter
using Plots
using Base.Meta: parse

# min |X / phi|^2 |PSD atm| + (|X / Lfast|^2 + |X / Lslow|^2) |PSD NCP| 
# subject to |Y / phi|^2 |PSD atm| + (|Y / Lfast|^2 + |Y / Lslow|^2) |PSD NCP| < threshold

f_cutoff = 15.0
f_loop = 1000.0
f_noise_crossover = 200.0
fr = 10 .^ (-2:0.01:log10(f_loop/2))
sr = 2π .* im .* fr
sT = sr / f_loop

ar1_high = ar1_filter(f_cutoff, f_loop, "high")
no_filter = ZPKFilter(0, 0, 1)
# f_loop, frame_delay, gain, leak, fpf

sys_high = AOSystem(f_loop, 1.0, 0.4, 0.999, 1, ar1_high)
sys_low = AOSystem(f_loop, 1.0, 1.4, 0.999, 1, no_filter)

Cfast = sT -> Hcont(sys_high, sT * f_loop) * Hfilter(sys_high, sT * f_loop)
Cslow = sT -> Hcont(sys_low, sT * f_loop)

vk_atm = VonKarman(v=10)
vk_ncp = VonKarman(v=0.1, rms_target=0.8)

begin
    plots = []
    for v in ["X", "Y"]
        ne = eval(parse("notched_error_$v"))
        p_v = plot(legend=:bottomright, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", ylims=(1e-10, 1e2), title="$v error = $(round(ne(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover), digits=3)) rad")
        for (fname, c, s) in zip(["phi_to_$v", "Lfast_to_$v", "Lslow_to_$v", "Nfast_to_$v", "Nslow_to_$v"], [1, 2, 2, 3, 3], [:solid, :solid, :dash, :solid, :dash])
            f = eval(parse(fname))
            plot!(fr, abs2.(f.(sT, Cfast, Cslow, 10)), label="|$fname|²", c=c, ls=s)
        end
        push!(plots, p_v)
    end
    plot(plots...)
end

