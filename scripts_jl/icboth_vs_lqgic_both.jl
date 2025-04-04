using multiwfs
using Plots

r0 = 0.1031
r0_ncp = 1.0
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))

f_loop = 1000.0
R = 10
leak = 0.999
no_filter = ZPKFilter(0, 0, 1)
fast_controller_ic = FilteredIntegrator(0.4, leak, ar1_filter(15.0, f_loop, "high"), 1/f_loop)
slow_controller_ic = FilteredIntegrator(1.4, leak, ZPKFilter(0, 0, 1), R/f_loop)
sim_ic = Simulation(f_loop, fast_controller_ic, slow_controller_ic, R, vk_atm, vk_ncp, 500.0)

nyquist_plot(sim_ic)

begin
    A = [[0. 0. 0. 0. 0. ]
    [1. 0. 0. 0. 0. ]
    [0. 0. leak 0. 0. ]
    [0. 0. leak 0. 0. ]
    [0. 0. 0. 1. 0. ]]
    B = [0.; 0; 1; 1; 0]
    C = [-1+1/R -1/R 0. 1-1/R 1/R]
    D = [1. 0 0 0 0]'
    Pw = hcat(1.0...)
    W = B * Pw * B'
    V = hcat(0.05...)
    K = kalman_gain(A, C, W, V)
    Vv = [0 1. 0 0 -1]
    Q = Vv' * Vv
    Rlqg = zeros(1,1)
    L = lqr_gain(A, D, Q, Rlqg)
    slow_controller_lqg = FilteredLQG(LQG(A, D, C, K, L, R/f_loop), no_filter)
end

begin
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
    fast_controller_lqg = FilteredLQG(LQG(Ã, D̃, C̃, K̃, L, 1/f_loop), ar1_filter(20.0, f_loop, "high"))
end

sim_lqgic = Simulation(f_loop, fast_controller_lqg, slow_controller_lqg, R, vk_atm, vk_ncp, 500.0)
nyquist_plot(sim_lqgic)

begin
    plot(sim_ic.fr, abs2.(phi_to_X.(sim_ic.sT, Ref(sim_ic))), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", label="IC, bandwidth = $(round(zero_db_bandwidth(sim_ic), digits=2)) Hz", legend=:bottomright, size=(400,400))
    plot!(sim_lqgic.fr, abs2.(phi_to_X.(sim_lqgic.sT, Ref(sim_lqgic))), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", label="LQG-IC, bandwidth = $(round(zero_db_bandwidth(sim_lqgic), digits=2)) Hz")
end

plot_integrands("X", sim_ic)
plot_integrands("X", sim_lqgic)