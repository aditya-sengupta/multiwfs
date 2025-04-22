using multiwfs
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * 0.6^(-5/3))

sim = simgen_ichpf(0.4, 0.0, 0.0, vk_ncp, 500.0);
plot_integrands("XY", sim)