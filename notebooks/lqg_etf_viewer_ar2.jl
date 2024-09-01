### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 6e63bc5c-6721-11ef-0e75-89cf19a9f817
begin
	using Pkg
	Pkg.activate("..")
	using multiwfs
	using Plots
	using LinearAlgebra: I
	using Base.GC: gc
	using PlutoUI
	using PlutoUI: combine
	using QuadGK
	using DataInterpolations
end

# ╔═╡ c57e0875-5d37-4693-9e6e-f02125cf4c0c
function design(params::Vector, lows::Vector, highs::Vector)
	
	return combine() do Child
		inputs = [
			md""" $(name): $(
				Child(name, Slider(l:0.001:h))
			)"""
			
			for (name, l, h) in zip(params, lows, highs)
		]
		
		md"""
		#### OL design
		$(inputs)
		"""
	end
end;

# ╔═╡ 98416ad3-5460-4b52-bf3b-073707a2e053
begin
	f_loop = 1000.0
	fr = exp10.(-4:0.01:log10(f_loop/2))
end;

# ╔═╡ 4622ed0c-6e04-4d18-a636-6316538c8b6c
vk = VonKarman();

# ╔═╡ fc73834e-37f1-4e3a-98e8-1d0e990fb3f2
@bind pars design(
	["freq_high", "damp_high", "freq_low", "damp_low", "log_hf_cost", "log_hf_process_noise"],
	[0.01, 0.00, 0.01, 0.00, -8.0, -8.0],
	[100.0, 1.0, 5.0, 1.0, 8.0, 8.0]
)

# ╔═╡ 0ae2ddcf-1341-412c-8560-7af718eb6e8b
function lqg_design_from_params(freq_high, damp_high, freq_low, damp_low, log_hf_cost, log_hf_process_noise, f_loop)
	Av1 = A_vib(freq_high/f_loop, damp_high)
    Av2 = A_vib(freq_low/f_loop, damp_low)
    L = A_DM(2)
    Ã = block_diag(L, Av1, Av2)
    C̃ = [0 -1 0 1 0 1]
    D̃ = [1 0 0 0 0 0]' 
    B = [0; 0; exp10(log_hf_process_noise); 0; 1; 0]
    Pw = hcat(1...)
    W = B * Pw * B'
    V = hcat(1...)
    K̃ = kalman_gain(Ã, C̃, W, V)
    Vv = [0 -1 0 exp10(log_hf_cost) 0 1]
    Q = Vv' * Vv
    R = zeros(1,1)
    L = lqr_gain(Ã, D̃, Q, R)
	return Ã, D̃, C̃, K̃, L
end

# ╔═╡ b84a3b69-c406-422c-a881-8b7858b6c5f3
function lqg_etf_from_params(freq_high, damp_high, freq_low, damp_low, log_hf_cost, log_hf_process_noise, f_loop)
	Ã, D̃, C̃, K̃, L = lqg_design_from_params(freq_high, damp_high, freq_low, damp_low, log_hf_cost, log_hf_process_noise, f_loop)
    s = 2π * im .* fr / f_loop
    z = exp.(s)
    zinvs = 1 ./ z
    gc()
	return 1 ./ (1 .+ lqg_controller_tf(Ã, D̃, C̃, K̃, L, zinvs))
end

# ╔═╡ 551dccd9-7348-4c18-8f9e-d2e65fba8a6f
begin
	lqg_etf_norm = abs2.(lqg_etf_from_params(pars.freq_high, pars.damp_high, pars.freq_low, pars.damp_low, pars.log_hf_cost, pars.log_hf_process_noise, f_loop))
	lqg_etf_norm_interp = CubicSpline(lqg_etf_norm, fr, extrapolate=true)
	residual_error = sqrt(quadgk(f -> psd_von_karman(f, vk) * lqg_etf_norm_interp(f), 0, 500)[1])
	plot(
		fr,
		lqg_etf_norm,
		xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="|ETF|²", xticks=[1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2],
		title="residual error = $(round(residual_error, digits=3)) rad", ylims=(1e-4, 1e4), label="ETF",
		legend=:topleft
	)
	vline!([pars.freq_low], label="freq_low = $(pars.freq_low) Hz, damp_low = $(pars.damp_low)", ls=:dash, color=:black)
	vline!([pars.freq_high], label="freq_high = $(pars.freq_high) Hz, damp_high = $(pars.damp_high)", ls=:dash, color=:black)
end

# ╔═╡ Cell order:
# ╠═6e63bc5c-6721-11ef-0e75-89cf19a9f817
# ╟─c57e0875-5d37-4693-9e6e-f02125cf4c0c
# ╟─98416ad3-5460-4b52-bf3b-073707a2e053
# ╟─4622ed0c-6e04-4d18-a636-6316538c8b6c
# ╟─fc73834e-37f1-4e3a-98e8-1d0e990fb3f2
# ╠═551dccd9-7348-4c18-8f9e-d2e65fba8a6f
# ╠═0ae2ddcf-1341-412c-8560-7af718eb6e8b
# ╠═b84a3b69-c406-422c-a881-8b7858b6c5f3
