### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 449bd402-0381-11f0-100f-19fe402916e9
begin
	using Pkg
	Pkg.activate("/Users/adityasengupta/projects/seal/multiwfs")
	using PlutoLinks
	using PlutoUI
	using Plots
	using Plots.PlotMeasures
end;

# ╔═╡ b418d4f6-a58d-4567-938b-702f9f755741
PlutoLinks.@revise using multiwfs

# ╔═╡ 4def56c6-90af-41c2-afbc-433b94d595b7
r0_ncp = 0.6

# ╔═╡ 609cfed9-0d6c-4a4f-b00e-339adab6459e
begin
	# atm/ncp setup
	r0 = 0.1031
	vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
	vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
	f_loop = 1000.0
	f_noise_crossover = 50.0
	R = 10
end;

# ╔═╡ 604e5a31-de48-4018-935f-a4eadb1b7440
function controllers_from_params(gain_slow, gain_fast, f_cutoff, leak=0.999)
	# controller setup
	fast_controller = FilteredIntegrator(gain_fast, leak, ar1_filter(f_cutoff, f_loop, "high"), 1/f_loop)
	slow_controller = FilteredIntegrator(gain_slow, leak, ZPKFilter(0, 0, 1), R/f_loop)
	return fast_controller, slow_controller
end;

# ╔═╡ 44e9e808-2bd2-4eaa-84f3-3c23f4ddb9be
begin
	fast_controller, slow_controller = controllers_from_params(1.4, 0.4, 15.0)
	sim = Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover)
end;

# ╔═╡ 1b5e2703-50b7-43e6-a519-5c730c0d0aa6
begin
	# hairdryer plot
	r0_ncp_vals = [0.6, 0.8, 1.0, 1.4, 1.6, 2.0, 4.0, 6.0]
	Xerrs_hdr = []
	Xerrs_hdr_nofilter = []
	fast_controller_nofilter = FilteredIntegrator(0.4, 0.999, ZPKFilter(0, 0, 1), 1/f_loop)
	slow_controller_off = FilteredIntegrator(0.0, 0.0, ZPKFilter(0, 0, 1), R/f_loop)
	for r0_ncp in r0_ncp_vals
		vk_ncp_hdr = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
		sim_hdr = Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp_hdr, f_noise_crossover)
		sim_hdr_nofilter = Simulation(f_loop, fast_controller_nofilter, slow_controller_off, R, vk_atm, vk_ncp_hdr, f_noise_crossover)
		search_gain!(sim_hdr_nofilter, "fast")
		push!(Xerrs_hdr, notched_error_X(sim_hdr))
		push!(Xerrs_hdr_nofilter, notched_error_X(sim_hdr_nofilter))
	end
	hairdryer = plot(r0_ncp_vals, Xerrs_hdr, xscale=:log10, xticks=(r0_ncp_vals, r0_ncp_vals), xlabel="NCP r₀ (m)", ylabel="X error (rad)", label="This controller", ylims=(0.5, 2.0))
	plot!(r0_ncp_vals, Xerrs_hdr_nofilter, xscale=:log10, xticks=(r0_ncp_vals, r0_ncp_vals), xlabel="NCP r₀ (m)", ylabel="X error (rad)", label="Integrator")
	vline!([r0_ncp], color=:black, ls=:dash, label="Reference r₀")
end;

# ╔═╡ 2d3b6440-3c01-483f-b9c4-b2f8a867817e
begin
	parameter_names = ["gain_slow", "gain_fast", "f_cutoff"]
	parameter_centers = [1.4, 0.4, 15.0]
	parameter_grids = [1.3:0.01:1.5, 0.3:0.01:0.5, 10.0:1.0:20.0]
	# resolved parameter minimum grid
	min_plots = []
	for (i, (parname, pargrid)) in enumerate(zip(parameter_names, parameter_grids))
		Xerrs_parsweep = []
		for par in pargrid
			fast_controller, slow_controller = controllers_from_params([i == j ? par : parameter_centers[j] for j in 1:3]...)
			sim_subopt = Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover)
			push!(Xerrs_parsweep, notched_error_X(sim_subopt))
		end
		p = plot(pargrid, Xerrs_parsweep, xlabel=parname, ylabel="X error", legend=nothing)
		vline!([parameter_centers[i]], color=:black, ls=:dash)
		push!(min_plots, p)
	end
end

# ╔═╡ 1837ea6d-9bfc-4976-afab-3481d07a98de
begin
	nyq = nyquist_plot(sim, legend=nothing)
	integrands = plot_integrands("XY", sim)
	etfplots = ten_etf_plots(sim)
	psdplot = five_psd_plots(sim)
end;

# ╔═╡ 206429bc-b933-4877-ae64-ad6107480839
begin
	allplots = []
	push!(allplots, hairdryer)
	push!(allplots, nyq)
	push!(allplots, integrands)
	append!(allplots, etfplots)
	push!(allplots, psdplot)
	append!(allplots, min_plots)
	pf = plot(allplots..., size=(1100,800), left_margin=5mm, suptitle="Lisa's IC-HPF; gain_slow=1.4, gain_fast=0.4, f_cutoff=15Hz, r0 NCP = $(r0_ncp)m", dpi=300)
	Plots.savefig(joinpath(multiwfs.PROJECT_ROOT, "figures", "evaluation", "evaluation_lisa_r0ncp$(r0_ncp).pdf"))
	pf
end

# ╔═╡ Cell order:
# ╠═449bd402-0381-11f0-100f-19fe402916e9
# ╠═b418d4f6-a58d-4567-938b-702f9f755741
# ╠═4def56c6-90af-41c2-afbc-433b94d595b7
# ╠═609cfed9-0d6c-4a4f-b00e-339adab6459e
# ╠═206429bc-b933-4877-ae64-ad6107480839
# ╠═604e5a31-de48-4018-935f-a4eadb1b7440
# ╠═44e9e808-2bd2-4eaa-84f3-3c23f4ddb9be
# ╠═1b5e2703-50b7-43e6-a519-5c730c0d0aa6
# ╟─2d3b6440-3c01-483f-b9c4-b2f8a867817e
# ╟─1837ea6d-9bfc-4976-afab-3481d07a98de
