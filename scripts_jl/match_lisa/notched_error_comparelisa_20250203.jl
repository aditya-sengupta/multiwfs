using multiwfs

f_loop = 1000.0
f_cutoff = 15.0
gain_slow, gain_fast = 1.4, 0.4
f_noise_crossover = 200.0
fr = 10 .^ (-2:0.01:log10(f_loop/2))
sr = 2π .* im .* fr
sT = sr / f_loop
ar1_high = ar1_filter(f_cutoff, f_loop, "high")
no_filter = ZPKFilter(0, 0, 1)
# f_loop, frame_delay, gain, leak, fpf
sys_slow = AOSystem(f_loop / 10, 0.1, gain_slow, 0.999, 1, ar1_high)
sys_fast = AOSystem(f_loop, 1.0, gain_fast, 0.999, 1, no_filter)
Cfast = sT -> Hcont(sys_fast, sT * f_loop) * Hfilter(sys_high, sT * f_loop)
Cslow = sT -> Hcont(sys_slow, sT * f_loop)
vk_atm = VonKarman(v=10)
vk_ncp = VonKarman(v=0.1, rms_target=0.8)
noise_normalization = psd_von_karman(200.0, vk_atm)

plot(fr, abs2.(Lfast_to_X.(sT, Cfast, Cslow, 10)), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", label="Lfast")
plot!(fr, abs2.(Nfast_to_X.(sT, Cfast, Cslow, 10)), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", label="Nfast", legend=:bottomright)

plot(fr, abs2.(Lslow_to_X.(sT, Cfast, Cslow, 10)), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", label="Lslow")
plot!(fr, abs2.(Nslow_to_X.(sT, Cfast, Cslow, 10)), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", label="Nslow", legend=:bottomright)

begin
    plots_etf = []
    for v in ["X", "Y"]
        ne = eval(parse("notched_error_$v"))
        p_v = plot(legend=:bottomright, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", ylims=(1e-10, 1e2), title="$v error = $(round(ne(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover), digits=3)) rad")
        for (fname, c, s) in zip(["phi_to_$v", "Lfast_to_$v", "Lslow_to_$v", "Nfast_to_$v", "Nslow_to_$v"], [1, 2, 2, 3, 3], [:solid, :solid, :dash, :solid, :dash])
            f = eval(parse(fname))
            plot!(fr, abs2.(f.(sT, Cfast, Cslow, 10)), label="|$fname|²", c=c, ls=s)
        end
        push!(plots_etf, p_v)
    end
    plot(plots_etf...)
end

begin
    plots = []
    for (v, sys) in zip(["X", "Y"], [sys_slow, [sys_slow, sys_fast]])
		tfc = is_stable(sys) ? :green : :red
        ne = eval(parse("notched_error_$v"))
        p_v = plot(legend=:bottomleft, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="Closed-loop residual at $v (rad)", title="\n $v error = $(round(ne(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover), digits=3)) rad \n", titlefontcolor=tfc, ylims=(1e-10, 1e2))
		for errsource in ["atm", "ncp", "noise"]
			err_source_fn = eval(parse("$(errsource)_error_at_f_$v"))
			plot!(fr, err_source_fn.(fr, Cfast, Cslow, 10, Ref(vk_atm), Ref(vk_ncp), noise_normalization, f_loop), label=errsource)
		end
		vline!([f_cutoff], color=:black, ls=:dash, label="HPF cutoff")
        push!(plots, p_v)
    end
    plot(plots..., suptitle="fc = $f_cutoff, gslow = $gain_slow, gfast = $gain_fast")
end