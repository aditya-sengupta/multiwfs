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

# ╔═╡ b4c74b15-516c-4faf-893c-5c211567d997
gain_slow = 1.52

# ╔═╡ 1a822cad-ba67-4e28-8446-b70a82d36e31
log_lqg_noise = 0.9

# ╔═╡ 4b8f47e2-829e-4c69-b403-5f7391b6c678
f_cutoff = 26.0

# ╔═╡ de44c458-82e3-45e9-8520-b8a588e923e2
num_nines_leak_slow = 1

# ╔═╡ 2d52af31-2c23-410e-a12c-1c469adce3b2
num_nines_leak_fast = 6

# ╔═╡ 82c0fbe2-1c4f-4fb4-8ebc-66779de6429b
leak_slow, leak_fast = 1 - exp10(-num_nines_leak_slow), 1 - exp10(-num_nines_leak_fast)

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
function controllers_from_params(gain_slow, log_lqg_noise, f_cutoff, leak_slow, leak_fast)
	A_ar1 = [leak_fast 0; 1 0]
	L = A_DM(2)
	Ã = block_diag(L, A_ar1)
	C̃ = [0 -1 0 1]
	D̃ = [1 0 0 0]' 
	B = [0; 0; 1; 0]
	Pw = hcat(1...)
	W = B * Pw * B'
	V = hcat(exp10(log_lqg_noise)...)
	K̃ = kalman_gain(Ã, C̃, W, V)
	Vv = [0 -1 0 1]
	Q = Vv' * Vv
	Rlqg = zeros(1,1)
	L = lqr_gain(Ã, D̃, Q, Rlqg)
	fast_controller = FilteredLQG(LQG(Ã, D̃, C̃, K̃, L, 1/f_loop), ar1_filter(f_cutoff, f_loop, "high"))
	slow_controller = FilteredIntegrator(gain_slow, leak_slow, ZPKFilter(0, 0, 1), R/f_loop)
	return fast_controller, slow_controller
end;

# ╔═╡ 44e9e808-2bd2-4eaa-84f3-3c23f4ddb9be
begin
	fast_controller, slow_controller = controllers_from_params(gain_slow, log_lqg_noise, f_cutoff, leak_slow, leak_fast)
	sim = Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover)
end;

# ╔═╡ 1b5e2703-50b7-43e6-a519-5c730c0d0aa6
begin
	# hairdryer plot
	r0_ncp_vals = [0.6, 0.8, 1.0, 1.4, 1.6, 2.0, 4.0, 6.0]
	Xerrs_hdr = []
	Xerrs_hdr_nofilter = []
	fast_controller_ref = FilteredIntegrator(0.4, 0.999, ZPKFilter(0, 0, 1), 1/f_loop)
	slow_controller_ref = FilteredIntegrator(1.4, 0.999, ZPKFilter(0, 0, 1), R/f_loop)
	for r0_ncp in r0_ncp_vals
		vk_ncp_hdr = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
		sim_hdr = Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp_hdr, f_noise_crossover)
		sim_hdr_nofilter = Simulation(f_loop, fast_controller_ref, slow_controller_ref, R, vk_atm, vk_ncp_hdr, f_noise_crossover)
		push!(Xerrs_hdr, notched_error_X(sim_hdr))
		push!(Xerrs_hdr_nofilter, notched_error_X(sim_hdr_nofilter))
	end
	hairdryer = plot(r0_ncp_vals, Xerrs_hdr, xscale=:log10, xticks=(r0_ncp_vals, r0_ncp_vals), xlabel="NCP r₀ (m)", ylabel="X error (rad)", label="This controller", ylims=(0.5, 2.0))
	plot!(r0_ncp_vals, Xerrs_hdr_nofilter, xscale=:log10, xticks=(r0_ncp_vals, r0_ncp_vals), xlabel="NCP r₀ (m)", ylabel="X error (rad)", label="(1.4, 0.4) integrator")
	vline!([r0_ncp], color=:black, ls=:dash, label="Reference r₀")
end

# ╔═╡ 055645ab-7f4e-4cf3-9945-f8391f6f485a
function nines_range(leak)
	nines = -log10(1 - leak)
	r = (nines-1):0.1:(nines+1)
	return @. 1 - exp10(-r)
end;

# ╔═╡ 2d3b6440-3c01-483f-b9c4-b2f8a867817e
begin
	parameter_names = ["gain_slow", "log_lqg_noise", "f_cutoff", "leak_slow", "leak_fast"]
	parameter_centers = [gain_slow, log_lqg_noise, f_cutoff, leak_slow, leak_fast]
	parameter_grids = [(gain_slow-0.1):0.01:(gain_slow+0.1), (log_lqg_noise-1):0.1:(log_lqg_noise+1), (f_cutoff-10):1:(f_cutoff+10), nines_range(leak_slow), nines_range(leak_fast)]
	# resolved parameter minimum grid
	min_plots = []
	for (i, (parname, pargrid)) in enumerate(zip(parameter_names, parameter_grids))
		Xerrs_parsweep = []
		for par in pargrid
			fast_controller, slow_controller = controllers_from_params([i == j ? par : parameter_centers[j] for j in 1:5]...)
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
	pf = plot(allplots..., size=(1100,1100), left_margin=5mm, suptitle="Leak-optimized fast-LQG-IC HPF; gain_slow=$gain_slow, fast LQG noise=$(round(exp10.(log_lqg_noise), digits=3)), f_cutoff=$f_cutoff \n leak_slow=$leak_slow, leak_fast=$leak_fast, r0 NCP = $(r0_ncp)m", dpi=300, layout=(4, 3))
	Plots.savefig(joinpath(multiwfs.PROJECT_ROOT, "figures", "evaluation", "evaluation_lqgicfast_hpf_leakopt_r0ncp$(r0_ncp).pdf"))
	pf
end

# ╔═╡ Cell order:
# ╠═449bd402-0381-11f0-100f-19fe402916e9
# ╠═b418d4f6-a58d-4567-938b-702f9f755741
# ╠═4def56c6-90af-41c2-afbc-433b94d595b7
# ╠═206429bc-b933-4877-ae64-ad6107480839
# ╠═b4c74b15-516c-4faf-893c-5c211567d997
# ╠═1a822cad-ba67-4e28-8446-b70a82d36e31
# ╠═4b8f47e2-829e-4c69-b403-5f7391b6c678
# ╠═de44c458-82e3-45e9-8520-b8a588e923e2
# ╠═2d52af31-2c23-410e-a12c-1c469adce3b2
# ╟─82c0fbe2-1c4f-4fb4-8ebc-66779de6429b
# ╠═609cfed9-0d6c-4a4f-b00e-339adab6459e
# ╠═604e5a31-de48-4018-935f-a4eadb1b7440
# ╠═44e9e808-2bd2-4eaa-84f3-3c23f4ddb9be
# ╠═1b5e2703-50b7-43e6-a519-5c730c0d0aa6
# ╠═055645ab-7f4e-4cf3-9945-f8391f6f485a
# ╠═2d3b6440-3c01-483f-b9c4-b2f8a867817e
# ╠═1837ea6d-9bfc-4976-afab-3481d07a98de
