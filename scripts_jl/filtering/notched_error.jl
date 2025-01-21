using multiwfs
using multiwfs: Hcont, Hfilter
using Plots
using Base.Meta: parse

# min |X / phi|^2 |PSD atm| + (|X / Lfast|^2 + |X / Lslow|^2) |PSD NCP| 
# subject to |Y / phi|^2 |PSD atm| + (|Y / Lfast|^2 + |Y / Lslow|^2) |PSD NCP| < threshold

f_cutoff = 2.0
f_loop = 1000.0
fr = 10 .^ (-2:0.01:log10(f_loop/20))
sr = 2π .* im .* fr / f_loop
zr = exp.(sr)

ar1_high = ar1_filter(f_cutoff, f_loop, "high")
no_filter = ZPKFilter(0, 0, 1)
# f_loop, frame_delay, gain, leak, fpf

sys_high = AOSystem(f_loop, 1.0, 0.4, 0.999, 10, ar1_high)
sys_low = AOSystem(f_loop / 10, 0.0, 1.7, 0.999, 1, no_filter)

Cfast = z -> Hcont(sys_high, log(z)) * Hfilter(sys_high, log(z))
Cslow = z -> Hcont(sys_low, log(z))

vk_atm = VonKarman(v=10)
vk_ncp = VonKarman(v=0.1, rms_target=0.8)

# first thing to look at, do my ETFs here match up with Lisa's?

notched_error_X(Cfast, Cslow, 10, vk_atm, vk_ncp)
notched_error_Y(Cfast, Cslow, 10, vk_atm, vk_ncp)

begin
    p_x = plot(legend=:left, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", ylims=(1e-10, 1e1), title="X error = $(round(notched_error_X(Cfast, Cslow, 10, vk_atm, vk_ncp), digits=3)) rad")
    for fname in ["phi_to_X", "Lfast_to_X", "Lslow_to_X"]
        f = eval(parse(fname))
        plot!(fr, abs2.(f.(zr, Cfast, Cslow, 10)), label="|$fname|²",)
    end
    p_y = plot(legend=:left, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", ylims=(1e-10, 1e1), title="Y error = $(round(notched_error_Y(Cfast, Cslow, 10, vk_atm, vk_ncp), digits=3)) rad")
    for fname in ["phi_to_Y", "Lfast_to_Y", "Lslow_to_Y"]
        f = eval(parse(fname))
        plot!(fr, abs2.(f.(zr, Cfast, Cslow, 10)), label="|$fname|²",)
    end
    plot(p_x, p_y)
end

plot(fr, abs2.(phi_to_X.(zr, Cfast, Cslow, 10)), xscale=:log10, yscale=:log10, label=nothing)

