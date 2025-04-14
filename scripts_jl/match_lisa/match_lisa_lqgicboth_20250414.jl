using multiwfs
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * 0.6^(-5/3))

begin
    simic = simgen_ichpf(0.4, 1.4, 15.0, vk_ncp, 500.0; leak=0.999)
    ten_etf_plots(simic, xlims=(0.01, 500.0), ylims=(1e-5, 1e1))[1]
end

begin
    sim = simgen_lqgicbothhpf(2.0, -1.0, 2.0, -1.0, 15.0, vk_ncp, 500.0);
    ten_etf_plots(sim, xlims=(0.01, 500.0), ylims=(1e-5, 1e1))[1]
end