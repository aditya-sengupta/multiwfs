### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 40749ad4-4447-4500-8a4d-ecd3ed8ea9ad
begin
	using Pkg
	Pkg.activate("..")
	using Plots
	using PlutoUI
	using multiwfs
end

# ╔═╡ 5d6a908f-209a-4555-af1e-8ec1331a70ab
num_nines_leak_ic = 1

# ╔═╡ f99012cc-7838-483a-b37d-d1ebad3bde54
@bind num_nines_leak_lqgic Slider(1:0.1:5)

# ╔═╡ de36d65b-d1e6-493e-9f35-c5ffa7eecf1d
begin
	leak_ic = 1 - exp10(-num_nines_leak_ic)
	leak_lqgic = 1 - exp10(-num_nines_leak_lqgic)
end;

# ╔═╡ f05718d3-7ec2-4f6a-a34a-aea718bef0ed
leak_ic

# ╔═╡ 85fb02f5-c424-4410-9cf5-4486b4899c81
begin
	f_loop = 1000.0
	R = 10
	vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * 0.1031^(-5/3))
	vk_ncp = VonKarman(0.3 * 0.01 / 100.0, 0.25)
	f_noise_crossover = 200.0
	no_filter = ZPKFilter(0, 0, 1)
	fast_controller = FilteredIntegrator(0.0, 0.0, no_filter, 1/f_loop)
	gain_ic = 1.655
	slow_controller = FilteredIntegrator(gain_ic, leak_ic, no_filter, R/f_loop)
	sim_ic = Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover)
end

# ╔═╡ 9e0f920f-4539-4930-9c84-4e4f2ffde6b3
begin
	delay_nframes = 1
	log_noise_lqgic = -8
    A = [[0. 0. 0. 0. 0. ]
    [1. 0. 0. 0. 0. ]
    [0. 0. leak_ic 0. 0. ]
    [0. 0. leak_ic 0. 0. ]
    [0. 0. 0. 1. 0. ]]
    B = [0.; 0; 1; 1; 0]
    C = [-(1-delay_nframes) -delay_nframes 0. 1-delay_nframes delay_nframes]
    D = [1. 0 0 0 0]'
    Pw = hcat(1.0...)
    W = B * Pw * B'
    V = hcat(exp10(log_noise_lqgic)...)
    K = kalman_gain(A, C, W, V)
    Vv = [0 1. 0 0 -1]
    Q = Vv' * Vv
    Rlqg = zeros(1,1)
    L = lqr_gain(A, D, Q, Rlqg)
    lqg = LQG(A, D, C, K, L, R/f_loop)
    slow_controller_lqg = FilteredLQG(lqg, no_filter)

    sim_lqgic = Simulation(f_loop, fast_controller, slow_controller_lqg, R, vk_atm, vk_ncp, f_noise_crossover)
end

# ╔═╡ d11e876d-2f75-4b84-8700-66535df5f174
begin
    plot(sim_ic.fr, abs2.(phi_to_X.(sim_ic.sT, Ref(sim_ic))), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", label="IC, bandwidth = $(round(zero_db_bandwidth(sim_ic), digits=2)) Hz", legend=:topleft, size=(400,400))
    plot!(sim_lqgic.fr, abs2.(phi_to_X.(sim_lqgic.sT, Ref(sim_lqgic))), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", label="LQG-IC, bandwidth = $(round(zero_db_bandwidth(sim_lqgic), digits=2)) Hz", legend=:topleft)
end

# ╔═╡ d148b5de-fb44-4eae-b465-18f8e70d7319
notched_error_X(sim_lqgic)

# add a noise component to the input PSD 
# and then find the optimal gain and leak for IC 
# And noise and alpha for LQG IC. 
# Adjust the noise level until the curve of gain vs. WFE is U-shaped 
# (i.e., the optimal gain is not due to stability limits). that case should clearly show a better LQG IC vs. IC.

# ╔═╡ 61323b57-8fa2-4e7a-9600-a72f45e94ed8
five_psd_plots(sim_ic)

# ╔═╡ 41688bc7-5b57-452b-bde7-03cf62b58e58


# ╔═╡ 5e01a3ff-2dbb-4d23-9b06-644e703499e0
begin
	notched_errors_ic = []
	gain_ic_evals = gain_ic-0.1:0.01:gain_ic+0.1
	for gain_ic_eval in gain_ic_evals
		slow_controller_eval = FilteredIntegrator(gain_ic_eval, leak_ic, no_filter, R/f_loop)
		sim_ic_eval = Simulation(f_loop, fast_controller, slow_controller_eval, R, vk_atm, vk_ncp, f_noise_crossover)
		push!(notched_errors_ic, notched_error_X(sim_ic_eval))
	end
	plot(gain_ic_evals, notched_errors_ic, xlabel="Gain (IC)", ylabel="X error", legend=nothing)
end

# ╔═╡ 828149d7-033f-4242-aabe-f5ce0e1b94e8
begin
	notched_errors_lqgic = []
	noise_lqgic_evals = exp10.(log_noise_lqgic-1:0.1:log_noise_lqgic+1)
	for noise_lqgic in noise_lqgic_evals
		A = [[0. 0. 0. 0. 0. ]
		    [1. 0. 0. 0. 0. ]
		    [0. 0. leak_ic 0. 0. ]
		    [0. 0. leak_ic 0. 0. ]
		    [0. 0. 0. 1. 0. ]]
		B = [0.; 0; 1; 1; 0]
		C = [-0.9 -0.1 0. 0.9 0.1]
		D = [1. 0 0 0 0]'
		Pw = hcat(1.0...)
		W = B * Pw * B'
		V = hcat(noise_lqgic...)
		K = kalman_gain(A, C, W, V)
		Vv = [0 1. 0 0 -1]
		Q = Vv' * Vv
		Rlqg = zeros(1,1)
		L = lqr_gain(A, D, Q, Rlqg)
		lqg = LQG(A, D, C, K, L, R/f_loop)
		slow_controller_lqg = FilteredLQG(lqg, no_filter)
		sim_lqgic_eval = Simulation(f_loop, fast_controller, slow_controller_lqg, R, vk_atm, vk_ncp, f_noise_crossover)
		push!(notched_errors_lqgic, notched_error_X(sim_lqgic_eval))
	end
	plot(noise_lqgic_evals, notched_errors_lqgic, xscale=:log10, xlabel="Noise (LQG-IC)", ylabel="X error", legend=nothing)
end


# ╔═╡ Cell order:
# ╠═40749ad4-4447-4500-8a4d-ecd3ed8ea9ad
# ╠═5d6a908f-209a-4555-af1e-8ec1331a70ab
# ╠═f05718d3-7ec2-4f6a-a34a-aea718bef0ed
# ╠═f99012cc-7838-483a-b37d-d1ebad3bde54
# ╟─d11e876d-2f75-4b84-8700-66535df5f174
# ╠═9e0f920f-4539-4930-9c84-4e4f2ffde6b3
# ╠═d148b5de-fb44-4eae-b465-18f8e70d7319
# ╟─de36d65b-d1e6-493e-9f35-c5ffa7eecf1d
# ╠═85fb02f5-c424-4410-9cf5-4486b4899c81
# ╠═61323b57-8fa2-4e7a-9600-a72f45e94ed8
# ╠═41688bc7-5b57-452b-bde7-03cf62b58e58
# ╠═5e01a3ff-2dbb-4d23-9b06-644e703499e0
# ╠═828149d7-033f-4242-aabe-f5ce0e1b94e8
