using multiwfs

r0 = 0.1031
r0_ncp = 1.0
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
f_loop = 1000.0
f_noise_crossover = 50.0
R = 10
leak = 0.999

function sim_generator_ic(gain_slow, gain_fast, f_cutoff)
    fast_controller = FilteredIntegrator(gain_fast, leak, ar1_filter(f_cutoff, f_loop, "high"), 1/f_loop)
	slow_controller = FilteredIntegrator(gain_slow, leak, no_filter, R/f_loop)
	return Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover)
end

grid_search_coarse_to_fine(
    sim_generator_ic,
    [[0.0, 2.0], [0.0, 1.0], [0.0, 100.0]];
    search="parallel"
)