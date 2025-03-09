using multiwfs
using Base.Meta: parse
using Plots

r0 = 0.1031
r0_ncp = 1.0
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))

# fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover
f_loop = 1000.0
R = 10
fast_controller = FilteredIntegrator(0.4, 0.999, ar1_filter(15.0, f_loop, "high"), 1/f_loop)
slow_controller = FilteredIntegrator(1.4, 0.999, ZPKFilter(0, 0, 1), R/f_loop)
sim = Simulation(f_loop, fast_controller, slow_controller, 10, vk_atm, vk_ncp, 500.0)

plot(sim.fr, abs2.(phi_to_X.(sim.sT, Ref(sim))), xscale=:log10, yscale=:log10)

notched_error_X(sim)
notched_error_Y(sim)

begin
    plots_etf = []
    for v in ["X", "Y"]
        ne = eval(parse("notched_error_$v"))
        p_v = plot(legend=:bottomright, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", ylims=(1e-10, 1e2), title="$v error = $(round(ne(sim), digits=3)) rad")
        for (fname, c, s) in zip(["Lslow_to_$v"], [1, 2, 2, 3, 3], [:solid, :solid, :dash, :solid, :dash])
            f = eval(parse(fname))
            plot!(sim.fr, abs2.(f.(sim.sT, Ref(sim))), label="|$fname|²", c=c, ls=s)
        end
        push!(plots_etf, p_v)
    end
    plot(plots_etf...)
end

begin
    errX, errY = notched_error_X(sim), notched_error_Y(sim)
    p_v = plot(legend=:bottomleft, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="Closed-loop residual (rad/Hz)", title="X error = $(round(errX, digits=3)) rad, Y error = $(round(errY, digits=3)) rad", ylims=(1e-10, 1e2))
    for errsource in ["atm_error_at_f_X", "ncp_error_at_f_X", "ncp_error_at_f_Y", "noise_error_at_f_X"]
        err_source_fn = eval(parse(errsource))
        plot!(fr, err_source_fn.(fr, Ref(sim)), label=errsource)
    end
    p_v
end