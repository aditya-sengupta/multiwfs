using multiwfs
using JLD2

# - pull an optimized IC-HPF and the equivalent LQG-IC-HPF

optpars = load("data/all_opt.jld2")
r0_ncp = 1.0
vk = VonKarman(0.001, 0.25 * r0_ncp^(-5/3))
f_crossover = 50.0

optpars["($r0_ncp, $f_crossover)"]

sim_ic = simgen_ichpf(optpars["($r0_ncp, $f_crossover)"][3]..., vk, f_crossover)
sim_lqgic = simgen_lqgicfasthpf(optpars["($r0_ncp, $f_crossover)"][4]..., vk, f_crossover)

notched_error_X(sim_ic)
notched_error_X(sim_lqgic)

# - take the LQG-IC parameters and compare just the $\phi$ to $X$ ETF for the fast controller between IC and LQG-IC. Is it better? Expect yes.

begin
    plot(sim_ic.fr, abs2.(phi_to_X.(sim_ic.sT, Ref(sim_ic))), xscale=:log10, yscale=:log10, label="Optimized IC", legend=:topleft)
    plot!(sim_lqgic.fr, abs2.(phi_to_X.(sim_lqgic.sT, Ref(sim_lqgic))), xscale=:log10, yscale=:log10, label="Optimized LQG-IC")
end

# - therefore, am I able to improve X error by using the f_cutoff from IC-HPF but the controller from LQG-IC? Expect something weird here.

sim_ic_without_cutoff = simgen_lqgicfasthpf(optpars["($r0_ncp, $f_crossover)"][3][1:end-1]..., 0.0, vk, f_crossover)

sim_lqgic_with_ichpfcutoff = simgen_lqgicfasthpf(optpars["($r0_ncp, $f_crossover)"][4][1:end-1]..., optpars["($r0_ncp, $f_crossover)"][3][3], vk, f_crossover)
notched_error_X(sim_lqgic_with_ichpfcutoff)

begin
    p1 = plot(sim_ic.fr, abs2.(phi_to_X.(sim_ic.sT, Ref(sim_ic))), xscale=:log10, yscale=:log10, label="Optimized IC", legend=:bottomright)
    plot!(sim_lqgic.fr, abs2.(phi_to_X.(sim_lqgic.sT, Ref(sim_lqgic))), xscale=:log10, yscale=:log10, label="Optimized LQG-IC")
    plot!(sim_lqgic_with_ichpfcutoff.fr, abs2.(phi_to_X.(sim_lqgic_with_ichpfcutoff.sT, Ref(sim_lqgic_with_ichpfcutoff))), xscale=:log10, yscale=:log10, label="LQG-IC plus filter", xlabel="Frequency (Hz)", ylabel="ETF atmosphere to X", xticks=exp10.(-3:2))
end

begin
    p2 = plot(sim_ic.fr, abs2.(Lslow_to_X.(sim_ic.sT, Ref(sim_ic))), xscale=:log10, yscale=:log10, label="Optimized IC", legend=:bottomright)
    plot!(sim_lqgic.fr, abs2.(Lslow_to_X.(sim_lqgic.sT, Ref(sim_lqgic))), xscale=:log10, yscale=:log10, label="Optimized LQG-IC")
    plot!(sim_lqgic_with_ichpfcutoff.fr, abs2.(Lslow_to_X.(sim_lqgic_with_ichpfcutoff.sT, Ref(sim_lqgic_with_ichpfcutoff))), xscale=:log10, yscale=:log10, label="LQG-IC plus filter", xlabel="Frequency (Hz)", ylabel="ETF X NCP to X", xticks=exp10.(-3:2))
end

plot(p1, p2,)