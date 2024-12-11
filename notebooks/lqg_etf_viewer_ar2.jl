### A Pluto.jl notebook ###
# v0.20.3

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

# ╔═╡ 6e63bc5c-6721-11ef-0e75-89cf19a9f817
begin
	using Pkg
	Pkg.activate("..")
	using multiwfs
	using multiwfs: Hfilter, Hcont, Lslow_to_Y
	using Plots
	using LinearAlgebra: I
	using Base.GC: gc
	using PlutoUI
	using PlutoUI: combine
	using QuadGK
	using DataInterpolations
	using Roots
	using Polynomials
end

# ╔═╡ 245f40df-5fb8-4152-9214-33a38cd21753
using LinearAlgebra: eigvals

# ╔═╡ c57e0875-5d37-4693-9e6e-f02125cf4c0c
function design(params::Vector, lows::Vector, highs::Vector)
	
	return combine() do Child
		inputs = [
			md""" $(name): $(
				Child(name, Slider(l:1e-4:h))
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
	fr = exp10.(-2:0.01:log10(f_loop/2))
	cutoff_idx = findfirst(x -> x > f_loop/20, fr)
	fr_slow = fr[1:cutoff_idx]
end;

# ╔═╡ 4622ed0c-6e04-4d18-a636-6316538c8b6c
vk = VonKarman();

# ╔═╡ 2cb452db-4559-49ee-9dcf-1ca164c433c8
"""begin
	freq_low = 0.0614
	damp_low = 1.0
	freq_high = 6.2523
	damp_high = 0.5666
	log_lf_cost = 0.2162
	log_lf_process_noise = -4.3804
	log_hf_cost = 1.001
	log_hf_process_noise = -5.3253
end;"""

# ╔═╡ fc73834e-37f1-4e3a-98e8-1d0e990fb3f2
@bind pars design(
	["ar1_param", "freq_low", "damp_low", "freq_high", "damp_high", "log_lf_cost", "log_lf_process_noise", "log_hf_cost", "log_hf_process_noise"],
	[0.5, 1e-4, 0.0, 1.0, 0.0, -8.0, -8.0, -8.0, -8.0],
	[0.999, 1e-1, 1.0, 100.0, 1.0, 8.0, 8.0, 8.0, 8.0]
)

# ╔═╡ 23f5e0d0-8f72-4791-b484-83f7a702b58a
ar1_param, freq_low, damp_low, freq_high, damp_high, log_lf_cost, log_lf_process_noise, log_hf_cost, log_hf_process_noise = pars.ar1_param, pars.freq_low, pars.damp_low, pars.freq_high, pars.damp_high, pars.log_lf_cost, pars.log_lf_process_noise, pars.log_hf_cost, pars.log_hf_process_noise

# ╔═╡ 9484c00f-2dde-4ddb-8b30-55947784025c
begin
	f_cutoff = 0.2
	ar1_low = ar1_filter(f_cutoff, f_loop / 10, "low")
	sys_low = AOSystem(f_loop, 1.0, 0.1, 0.9999999, 10, ar1_low)
	search_gain!(sys_low)
	lpf_ol = Hol.(Ref(sys_low), fr_slow)
	low_etf = Hrej.(Ref(sys_low), fr_slow)
end;

# ╔═╡ 0ae2ddcf-1341-412c-8560-7af718eb6e8b
function lqg_design_from_params(ar1_param, freq_high, damp_high, freq_low, damp_low, log_lf_cost, log_lf_process_noise, log_hf_cost, log_hf_process_noise, f_loop)
	Av1 = A_vib(freq_high/f_loop, damp_high)
    Av2 = A_vib(freq_low/f_loop, damp_low)
	A_ar1 = [ar1_param 0; 1 0]
    L = A_DM(2)
    Ã = block_diag(L, A_ar1, Av1)
    C̃ = [0 -1 0 1 0 1]
    D̃ = [1 0 0 0 0 0]' 
    B = [0; 0; 1; 0; exp10(log_hf_process_noise); 0]
    Pw = hcat(1...)
    W = B * Pw * B'
    V = hcat(1...)
    K̃ = kalman_gain(Ã, C̃, W, V)
    Vv = [0 -1 0 1 0 exp10(log_hf_cost)]
    Q = Vv' * Vv
    R = zeros(1,1)
    L = lqr_gain(Ã, D̃, Q, R)
	return Ã, D̃, C̃, K̃, L
end;

# ╔═╡ a12e3c9f-efd1-4ac2-9825-2f8aeb3afb7f
lqg = LQG(lqg_design_from_params(ar1_param, freq_high, damp_high, freq_low, damp_low, log_lf_cost, log_lf_process_noise, log_hf_cost, log_hf_process_noise, f_loop)...);

# ╔═╡ feaa916d-d741-464f-807a-cfca2f8cc208
begin
	sr = 2π * im * fr ./ f_loop
	zr = exp.(sr)
	Cfast = z -> transfer_function.(Ref(lqg), log(z))
	Cslow = z -> (real(log(z) / (2π * im)) < 1/20) ? (Hfilter(sys_low, log(z)) * Hcont(sys_low, log(z))) : 1.0
end;

# ╔═╡ 7896c60d-5fc4-410f-9f4a-816ed611c149
eigvals(lqg.ikcA) .|> abs |> maximum

# ╔═╡ b84a3b69-c406-422c-a881-8b7858b6c5f3
function lqg_etf_from_params(ar1_param, freq_high, damp_high, freq_low, damp_low, log_lf_cost, log_lf_process_noise, log_hf_cost, log_hf_process_noise, f_loop)
	Ã, D̃, C̃, K̃, L = lqg_design_from_params(ar1_param, freq_high, damp_high, freq_low, damp_low, log_lf_cost, log_lf_process_noise, log_hf_cost, log_hf_process_noise, f_loop)
    s = 2π * im .* fr / f_loop
    z = exp.(s)
    zinvs = 1 ./ z
    gc()
	return 1 ./ (1 .+ lqg_controller_tf(Ã, D̃, C̃, K̃, L, zinvs))
end;

# ╔═╡ 551dccd9-7348-4c18-8f9e-d2e65fba8a6f
begin
	lqg_etf = lqg_etf_from_params(ar1_param, freq_high, damp_high, freq_low, damp_low, log_lf_cost, log_lf_process_noise, log_hf_cost, log_hf_process_noise, f_loop)
	lqg_etf_norm = abs2.(lqg_etf)
	lqg_etf_norm_interp = CubicSpline(lqg_etf_norm, fr, extrapolate=true)
	residual_error = sqrt(quadgk(f -> psd_von_karman(f, vk) * lqg_etf_norm_interp(f), 0.1, 500)[1])
	hpf_ol = 1 ./ lqg_etf .- 1
	p = plot(
		fr,
		lqg_etf_norm,
		xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="|ETF|²", xticks=[1e-2, 1e-1, 1e0, 1e1, 1e2], ylims=(1e-8, 1e2), label="High-speed ETF",
		legend=:bottomright
	)
	plot!(
		fr_slow,
		abs2.(low_etf),
		xscale=:log10, yscale=:log10, label="Low-speed ETF",
	)
	combined_etf_norm = abs2.(1 ./ (1 .+ hpf_ol[1:cutoff_idx] .+ lpf_ol))
	append!(combined_etf_norm, lqg_etf_norm[cutoff_idx+1:end])
	plot!(
		fr,
		combined_etf_norm,
		lw=2, color=4, label="Combined ETF"
	)
	hline!([1], ls=:dash, color=:black, label=nothing)
	p2 = plot(fr, Lfast_to_X.(zr, Cfast, Cslow, 10) .|> abs2, label="Fast NCP -> Slow WFS", xscale=:log10, yscale=:log10, ylims=(1e-8, 1e2), legend=:bottomright, xlabel="Frequency (Hz)")
	plot!(fr, Lslow_to_Y.(zr, Cfast, Cslow, 10) .|> abs2, label="Slow NCP -> Fast WFS")
	#plot!(fr, Y_to_X.(zr, Cfast, Cslow, 10) .|> abs2, color=:green, label="Fast WFS -> Slow WFS")
	hline!([1], ls=:dash, color=:black, label=nothing)
	# Plots.savefig("../figures/lqgfirst_etf_parts.pdf")
	plot(p, p2, suptitle="HPF res = $(round(residual_error, digits=3)) rad, max A eig = $(round(maximum(abs.(eigvals(lqg.A))), digits=5))")
end

# ╔═╡ 9528297f-875e-4e7f-9592-61bcac0dd197
"""npzwrite(
    "../data/lqgfirst_limitmaxeig.npz",
    Dict(
        "K" => lqg.K,
        "G" => lqg.L,
        "ikcA" => lqg.ikcA,
        "ikcD" => lqg.ikcD
    )
)"""

# ╔═╡ 4cd07eca-976c-4230-8031-45a2004bbbbd
begin
	A_low_AR = [0.995 0; 1 0]
	L = A_DM(2)
	A_low = block_diag(L, A_low_AR)
	D_low = [1 0 0 0]'
	C_low = [0 -1 0 1]
	B_low = [0; 0; 1; 0]
	W_low = B_low * hcat(1...) * B_low'
	V_low = hcat(1...)
	K_low = kalman_gain(A_low, C_low, W_low, V_low)
	Vv_low = [0 -1 0 1]
	Q_low = Vv_low' * Vv_low
	R_low = zeros(1,1)
	L_low = lqr_gain(A_low, D_low, Q_low, R_low)
	s = 2π * im .* fr_slow / (f_loop / 10)
    z = exp.(s)
    zinvs = 1 ./ z
	L
	low_etf_lqg = 1 ./ (1 .+ lqg_controller_tf(A_low, D_low, C_low, K_low, L_low, zinvs))
end;

# ╔═╡ Cell order:
# ╠═6e63bc5c-6721-11ef-0e75-89cf19a9f817
# ╟─c57e0875-5d37-4693-9e6e-f02125cf4c0c
# ╟─98416ad3-5460-4b52-bf3b-073707a2e053
# ╠═4622ed0c-6e04-4d18-a636-6316538c8b6c
# ╠═2cb452db-4559-49ee-9dcf-1ca164c433c8
# ╠═23f5e0d0-8f72-4791-b484-83f7a702b58a
# ╠═fc73834e-37f1-4e3a-98e8-1d0e990fb3f2
# ╟─feaa916d-d741-464f-807a-cfca2f8cc208
# ╠═551dccd9-7348-4c18-8f9e-d2e65fba8a6f
# ╠═245f40df-5fb8-4152-9214-33a38cd21753
# ╠═7896c60d-5fc4-410f-9f4a-816ed611c149
# ╠═9484c00f-2dde-4ddb-8b30-55947784025c
# ╠═0ae2ddcf-1341-412c-8560-7af718eb6e8b
# ╟─a12e3c9f-efd1-4ac2-9825-2f8aeb3afb7f
# ╠═b84a3b69-c406-422c-a881-8b7858b6c5f3
# ╠═9528297f-875e-4e7f-9592-61bcac0dd197
# ╟─4cd07eca-976c-4230-8031-45a2004bbbbd
