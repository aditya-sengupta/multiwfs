using multiwfs

r0 = 0.1031
r0_ncp = 0.6
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
f_loop = 1000.0
f_noise_crossover = 50.0
R = 10

leak = 0.995

    fast_controller = FilteredIntegrator(0.4, leak, ar1_filter(15.0, f_loop, "high"), 1/f_loop)
    slow_controller = FilteredIntegrator(1.4, leak, ZPKFilter(0, 0, 1), R/f_loop)
    sim = Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover)

plot(ten_etf_plots(sim)[1], ylim=(1e-5,1e1), xlim=(0.0125, 500.0))