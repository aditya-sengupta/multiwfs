using multiwfs
using JLD2
using Plots
using Colors

f_crossover = 200.0

begin
    optpars = load("data/all_opt.jld2")
    xerrs = Dict()
    r0_ncps = vcat(0.6:0.2:2.0, 3:6)
    for r0_ncp in r0_ncps
        vk_ncp = VonKarman(0.001, 0.25 * (r0_ncp)^(-5/3))
        sim_one = simgen_ichpf(optpars["($r0_ncp, $f_crossover)"][1][1], 0.0, 0.0, vk_ncp, f_crossover; leak=0.999)
        sim_two = simgen_ichpf(optpars["($r0_ncp, $f_crossover)"][3]..., vk_ncp, f_crossover; leak=0.999)
        xerrs["($r0_ncp, $f_crossover) one"] = notched_error_X(sim_one)
        xerrs["($r0_ncp, $f_crossover) one noncp"] = notched_error_X_noncp(sim_one)
        xerrs["($r0_ncp, $f_crossover) two"] = notched_error_X(sim_two)
        xerrs["($r0_ncp, $f_crossover) two noncp"] = notched_error_X_noncp(sim_two)
    end
    
    p = Plots.plot(xlabel="r₀ NCP (m)", ylabel="CL error at X (rad)", xscale=:log10, xticks=([0.6, 1.0, 1.4, 2.0, 3.0, 4.0, 6.0], [0.6, 1.0, 1.4, 2.0, 3.0, 4.0, 6.0]))
    Plots.plot!(r0_ncps, [xerrs["($r0_ncp, $f_crossover) one noncp"] for r0_ncp in r0_ncps], msw=0, label="Main: Atm only", color=RGBA(1/255,115/255,178/255,255/255), lw=2)
    Plots.plot!(r0_ncps, [xerrs["($r0_ncp, $f_crossover) two noncp"] for r0_ncp in r0_ncps], msw=0, label="Both: Atm only", color=RGBA(222/255,143/255,5/255,255/255), lw=2)
    Plots.plot!(r0_ncps, [xerrs["($r0_ncp, $f_crossover) one"] for r0_ncp in r0_ncps], msw=0, label="Main: Atm + NCP", color=RGBA(2/255,158/255,115/255,255/255), lw=2)
    Plots.plot!(r0_ncps, [xerrs["($r0_ncp, $f_crossover) two"] for r0_ncp in r0_ncps], msw=0, label="Both: Atm + NCP", color=RGBA(213/255,94/255,0/255,255/255), lw=2)
    Plots.plot(p)
end

# r0 3m: check that the control p
threem_pars = optpars["(3.0, 200.0)"]
threem_pars[1][1] 
# fast gain for one integrator
threem_pars[3]
# fast gain, slow gain, cutoff freq for two integrators + HPF

vk_ncp = VonKarman(0.001, 0.25 * 3^(-5/3))
sim_one = simgen_ichpf(threem_pars[1][1], 0.0, 0.0, vk_ncp, f_crossover; leak=0.999)
sim_two = simgen_ichpf(threem_pars[3]..., vk_ncp, f_crossover; leak=0.999)

threem_pars[3]

ten_etf_plots(sim_one, xlim=(1,500), ylim=(1e-4,1e2), suptitle="one WFS ETFs")[1]
ten_etf_plots(sim_two, xlim=(1,500), ylim=(1e-4,1e2), suptitle="two WFS ETFs")[1]

plot_integrands("XY", sim_one, 1:4, xlim=(1,500), suptitle="one WFS integrands")
plot_integrands("XY", sim_two, 1:4, xlim=(1,500), suptitle="two WFS integrands")

zero_db_bandwidth(sim_one)
zero_db_bandwidth(sim_two)

atm_error_at_f_X(100.0, sim_one)
atm_error_at_f_X(100.0, sim_two)

ncp_error_at_f_X(100.0, sim_one)
ncp_error_at_f_X(100.0, sim_two)

noise_error_at_f_X(100.0, sim_one)
noise_error_at_f_X(100.0, sim_two)