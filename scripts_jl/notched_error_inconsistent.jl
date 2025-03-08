using multiwfs

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

notched_error_X(sim)
notched_error_Y(sim)