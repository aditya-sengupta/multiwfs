### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 449bd402-0381-11f0-100f-19fe402916e9
begin
	using Pkg
	Pkg.activate("../..")
	using PlutoLinks
	using Plots
end

# ╔═╡ b418d4f6-a58d-4567-938b-702f9f755741
PlutoLinks.@revise using multiwfs

# ╔═╡ 604e5a31-de48-4018-935f-a4eadb1b7440
begin
	# simulation setup
	design = "lisa"
	r0 = 0.1031
	r0_ncp = 0.6
	vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
	vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
	f_loop = 1000.0
	f_noise_crossover = 500.0
	R = 10
	leak = 0.999
	fast_controller = FilteredIntegrator(0.4, leak, ar1_filter(15.0, f_loop, "high"), 1/f_loop)
	slow_controller = FilteredIntegrator(1.4, leak, ZPKFilter(0, 0, 1), R/f_loop)
	sim = Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover)
end;

# ╔═╡ 1b5e2703-50b7-43e6-a519-5c730c0d0aa6
begin
	# hairdryer plot
	r0_ncp_vals = [0.6, 0.8, 1.0, 1.4, 1.6, 2.0, 4.0, 6.0]
	Xerrs_hdr = []
	for r0_ncp in r0_ncp_vals
		vk_ncp_hdr = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
		sim_hdr = Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp_hdr, f_noise_crossover)
		push!(Xerrs_hdr, notched_error_X(sim_hdr))
	end
	plot(r0_ncp_vals, Xerrs_hdr, xscale=:log10, xticks=(r0_ncp_vals, r0_ncp_vals), xlabel="NCP r₀ (m)", ylabel="X error (rad)", label=nothing, ylims=(0.5, 1.0), title="Hairdryer plot, $design")
end

# ╔═╡ 24f7b846-3ec1-480d-af07-becda7dd8ee0
nyquist_plot(sim, title="Nyquist plot, $design")

# ╔═╡ 66674200-839a-4816-a169-60dae3ab7d3e
plot_integrands(sim, legend=:topleft)

# ╔═╡ Cell order:
# ╠═449bd402-0381-11f0-100f-19fe402916e9
# ╠═b418d4f6-a58d-4567-938b-702f9f755741
# ╟─604e5a31-de48-4018-935f-a4eadb1b7440
# ╟─1b5e2703-50b7-43e6-a519-5c730c0d0aa6
# ╟─24f7b846-3ec1-480d-af07-becda7dd8ee0
# ╠═66674200-839a-4816-a169-60dae3ab7d3e
