### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# ╔═╡ 8f38fa65-b607-4ca2-bcc6-1afb6eb144ce
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 60ba798b-7753-4080-8c6d-0c31e069ac45
begin
	using CairoMakie
	using StatsBase: mean
	using PlutoLinks
	using JLD2
	using NPZ
	PlutoLinks.@revise using multiwfs
end

# ╔═╡ c2a9ade5-9c78-4f02-91ba-f8e327c9ebab
begin
	optpars_ichpf = load("../data/all_opt.jld2")
	optpars_ichpf = Dict(k => optpars_ichpf[k][3] for k in keys(optpars_ichpf))
end 

# ╔═╡ da2ec842-5c06-4bdb-8d2b-64438de1e010
r0_ncp = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, 6.0]

# ╔═╡ 95a38ad5-ac79-47a1-8916-b4d36b8ca7e9
optpars_ichpf["(1.0, 50.0)"]

# ╔═╡ 5c2dac77-f2d8-45d1-b6eb-34dbb62200b6
begin
	#olt = generate_openloop_timeseries(sim_ref, 100_000)
	#npzwrite("../data/olt/olt_atm.npy", olt.atm)
	#npzwrite("../data/olt/olt_fastncp.npy", olt.fastncp)
	#npzwrite("../data/olt/olt_slowncp.npy", olt.slowncp)
	#npzwrite("../data/olt/olt_fastnoise.npy", olt.fastnoise)
	#npzwrite("../data/olt/olt_slownoise.npy", olt.slownoise)
end

# ╔═╡ ec9b4b81-61a7-476a-bd08-d3d9b1892234
begin
	N = 100_000
	vk_ncp = VonKarman(0.001, 0.25 * (1)^(-5/3))
	f_crossover = 500.0
	sim = simgen_ichpf(0.4, 1.4, 15.0, vk_ncp, f_crossover; leak=0.995)
	all_timeseries = generate_openloop_timeseries(sim, N)
end;

# ╔═╡ 8b2008be-e6b7-4518-b4ea-aae71f5ab20c
sim_ref = simgen_ichpf(optpars_ichpf["(1.0, 50.0)"]..., VonKarman(0.001, 0.25), f_crossover)

# ╔═╡ 9f3111cb-2de1-4ea4-812b-ce3aed6f19bc
begin
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="r₀ NCP (m)", ylabel="CL error at X (rad)", xscale=log10, xticks=r0_ncp)
	for (f_crossover, c) in zip([50.0, 100.0, 200.0, 500.0], [:blue, :green, :red, :purple])
		rms_at_X = []
		notched_error_freqdomain = []
		for r0_ncp_val in r0_ncp
			local vk_ncp = VonKarman(0.001, 0.25 * r0_ncp_val^(-5/3))
			sim = simgen_ichpf(optpars_ichpf["($r0_ncp_val, $f_crossover)"]..., vk_ncp, f_crossover)
			this_rms_at_X = 0.0
			Ncases = 1
			for _ in 1:Ncases
				all_timeseries = generate_openloop_timeseries(sim, N)
				closedloop_at_X = timedomain_closedloop(sim, all_timeseries)
				this_rms_at_X += rms(closedloop_at_X) / Ncases
			end
			push!(rms_at_X, this_rms_at_X)
			push!(notched_error_freqdomain, notched_error_X(sim))
		end
		lines!(ax, r0_ncp, rms_at_X, label="f crossover = $f_crossover Hz", color=c)
		lines!(ax, r0_ncp, notched_error_freqdomain, color=c, linestyle=:dash)
	end
	axislegend(ax)
	fig
end

# ╔═╡ 75bf5a1e-8468-4972-b4ff-0ac4b11a0e8b
closedloop_at_X = timedomain_closedloop(sim, all_timeseries);

# ╔═╡ f300cdcc-51b7-4403-98de-f0b3f9f55769
plot_timedomain_vs_analytic(sim, phi_to_X, all_timeseries.atm, closedloop_at_X)

# ╔═╡ e2c676b3-23dc-4793-98dd-224a583bc32d
plot_timedomain_vs_analytic(sim, Lfast_to_X, all_timeseries.fastncp, closedloop_at_X)

# ╔═╡ d11f52a7-157b-40d6-89df-e086b3b3a0d2
plot_timedomain_vs_analytic(sim, Lslow_to_X, all_timeseries.slowncp, closedloop_at_X)

# ╔═╡ 182a422d-9f57-4988-a7a1-9379b914d76a
plot_timedomain_vs_analytic(sim, Nfast_to_X, all_timeseries.fastnoise, closedloop_at_X)

# ╔═╡ 49f2608a-9bce-416c-9468-aa865ed6f965
plot_timedomain_vs_analytic(sim, Nslow_to_X, all_timeseries.slownoise, closedloop_at_X)

# ╔═╡ Cell order:
# ╠═8f38fa65-b607-4ca2-bcc6-1afb6eb144ce
# ╠═60ba798b-7753-4080-8c6d-0c31e069ac45
# ╠═c2a9ade5-9c78-4f02-91ba-f8e327c9ebab
# ╠═da2ec842-5c06-4bdb-8d2b-64438de1e010
# ╠═95a38ad5-ac79-47a1-8916-b4d36b8ca7e9
# ╠═8b2008be-e6b7-4518-b4ea-aae71f5ab20c
# ╠═5c2dac77-f2d8-45d1-b6eb-34dbb62200b6
# ╠═9f3111cb-2de1-4ea4-812b-ce3aed6f19bc
# ╠═ec9b4b81-61a7-476a-bd08-d3d9b1892234
# ╠═75bf5a1e-8468-4972-b4ff-0ac4b11a0e8b
# ╠═f300cdcc-51b7-4403-98de-f0b3f9f55769
# ╠═e2c676b3-23dc-4793-98dd-224a583bc32d
# ╠═d11f52a7-157b-40d6-89df-e086b3b3a0d2
# ╠═182a422d-9f57-4988-a7a1-9379b914d76a
# ╠═49f2608a-9bce-416c-9468-aa865ed6f965
