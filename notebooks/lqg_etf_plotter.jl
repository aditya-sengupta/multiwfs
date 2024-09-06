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

# ╔═╡ d21f40e9-24ac-4830-9f82-7a23f6994432
begin
	using Pkg
	Pkg.activate(".")
	using LinearAlgebra: I
	using Plots
	using multiwfs: A_vib, block_diag, kalman_gain, lqr_gain
	using PlutoUI
	using PlutoUI: combine
end;

# ╔═╡ 3fc2059a-2785-4ee0-b186-b8c918f474e0
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
end

# ╔═╡ ecf4806d-a3f5-4f31-b5eb-ebbec17db0a8
@bind pars design(
	["freq_high", "damp_high", "freq_low", "damp_low", "log_hf_cost"],
	[10.0, 0.00, 0.01, 0.00, -7.0],
	[500.0, 1.0, 10.0, 1.0, 7.0]
)

# ╔═╡ 46acb1d2-347d-4686-a3a6-94f513fc15a4
begin
	f_loop = 1000.0
	fr = exp10.(-4:0.01:log10(f_loop/2))
	s = 2π * im * fr ./ f_loop
	z = exp.(s)
	zinvs = 1 ./ z
	n_input_history = 3
	L = zeros(n_input_history, n_input_history)
	for i in 1:(n_input_history-1)
		L[i+1,i] = 1
	end
	C = [0 0 -1 1 0 1 0]
	B̃ = [1 0 0 0 0 0 0]'
	V = hcat(8.0...)
	W_g = zeros(7,7)
	R = zeros(1,1)
end;

# ╔═╡ 58b7b95c-4114-435f-ae1b-c4dff41a3b8a
function lqg_tf(A, D, C, K, G)
    ikcA = (I - K * C) * A
    ikcD = (I - K * C) * D
    numerator = [(G * inv(I - ikcA * zinv))[1,1] for zinv in zinvs]
	denominator = [(I - G * inv(I - ikcA * zinv) * ikcD * zinv)[1,1] for zinv in zinvs]
	plant = numerator .* (zinvs .^ 2) ./ denominator
    return 1 ./ (1 .+ plant)
end;

# ╔═╡ a722afcd-0204-4dd2-a8b0-9f11eb0f9e3d
function lqg_tf_from_setup(freq_high, damp_high, freq_low, damp_low, highfreq_cost, w)
    Av1 = real.(A_vib(freq_high/f_loop, damp_high))
    Av2 = real.(A_vib(freq_low/f_loop, damp_low))
	A = block_diag(L, Av1, Av2)
	Ccost = [1 0 0 highfreq_cost 0 1 0]
	W = copy(W_g)
	W[n_input_history+1:end,n_input_history+1:end] = w * I(4)
	Q = Ccost' * Ccost
	K = kalman_gain(A, C, W, V)
	G = lqr_gain(A, B̃, Q, R)
	lqgf = lqg_tf(A, B̃, C, K, G)
end;

# ╔═╡ a7c24efb-0d32-4103-864c-9a38769cd883
highfreq_cost = exp10(pars.log_hf_cost);

# ╔═╡ acb46173-dbdb-4c47-8cf5-5c58c68a7231
begin
	pl = plot(legend=:topleft)
	ws = exp10.(-0:2:0)
	for (i, w) in enumerate(ws)
		lqgf = nothing
		try
			lqgf = lqg_tf_from_setup(pars.freq_high, pars.damp_high, pars.freq_low, pars.damp_low, highfreq_cost, w)
		catch e
			lqgf = lqg_tf_from_setup(pars.freq_high, pars.damp_high + 0.01, pars.freq_low, pars.damp_low, highfreq_cost, w)
		end
		plot!(fr, abs2.(lqgf), xscale=:log10, yscale=:log10, label="w = $w", xlabel="Frequency (Hz)", ylabel="|ETF|²", title="freqs = ($(pars.freq_low), $(pars.freq_high)), damp = ($(pars.damp_low), $(pars.damp_high)), log10 HF cost = $(pars.log_hf_cost)", legend=:bottomleft, titlefontsize=10, xticks=[1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2], ylims=(1e-3, 10))
	end
	vline!([f_loop/2], label="Nyquist limit", ls=:dash, color=:black)
	hline!([1], label="Unit rejection", ls=:dash, color=:red)
	pl
end

# ╔═╡ a6694103-e5d7-49bf-87c4-84c388eef396
lqgf = lqg_tf_from_setup(pars.freq_high, pars.damp_high, pars.freq_low, pars.damp_low, highfreq_cost, 1)

# ╔═╡ Cell order:
# ╠═d21f40e9-24ac-4830-9f82-7a23f6994432
# ╠═58b7b95c-4114-435f-ae1b-c4dff41a3b8a
# ╟─3fc2059a-2785-4ee0-b186-b8c918f474e0
# ╟─ecf4806d-a3f5-4f31-b5eb-ebbec17db0a8
# ╠═acb46173-dbdb-4c47-8cf5-5c58c68a7231
# ╠═46acb1d2-347d-4686-a3a6-94f513fc15a4
# ╠═a722afcd-0204-4dd2-a8b0-9f11eb0f9e3d
# ╟─a6694103-e5d7-49bf-87c4-84c388eef396
# ╟─a7c24efb-0d32-4103-864c-9a38769cd883
