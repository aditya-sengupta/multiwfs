using multiwfs
using Plots

r0 = 0.1031
r0_ncp = 1.0
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))

# fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover
f_loop = 1000.0
R = 10
fast_controller = FilteredIntegrator(0.43, 0.999, no_filter, 1/f_loop)
slow_controller = FilteredIntegrator(0.0, 0.999, no_filter, R/f_loop)
sim = Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, 500.0)

A_ar1 = [0.999 0; 1 0]
L = A_DM(2)
Ã = block_diag(L, A_ar1)
C̃ = [0 -1 0 1]
D̃ = [1 0 0 0]' 
B = [0; 0; 1; 0]
Pw = hcat(1...)
W = B * Pw * B'
V = hcat(0.05...)
K̃ = kalman_gain(Ã, C̃, W, V)
Vv = [0 -1 0 1]
Q = Vv' * Vv
Rlqg = zeros(1,1)
L = lqr_gain(Ã, D̃, Q, Rlqg)
lqgic_controller = LQG(Ã, D̃, C̃, K̃, L, 1/sim.f_loop)
sim_lqgic = Simulation(f_loop, lqgic_controller, slow_controller, R, vk_atm, vk_ncp, 500.0)

lqgic_tf_delay = abs2.(phi_to_X.(sim.sT, Ref(sim_lqgic)))
bw_delay = zero_db_bandwidth(sim_lqgic)

begin
    plot(sim.fr, abs2.(phi_to_X.(sim.sT, Ref(sim))), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", label="IC, bandwidth = $(round(zero_db_bandwidth(sim), digits=2)) Hz", legend=:bottomleft, size=(400,400))
    plot!(sim.fr, abs2.(phi_to_X.(sim.sT, Ref(sim_lqgic))), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", label="LQG-IC, bandwidth = $(round(zero_db_bandwidth(sim_lqgic), digits=2)) Hz")
end
