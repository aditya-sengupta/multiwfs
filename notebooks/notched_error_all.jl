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

# ╔═╡ 65e6fa3f-0342-434b-b37d-8a0ec7ddc1f9
begin
	using Pkg
	Pkg.activate("..")
	using multiwfs
end

# ╔═╡ 247e2d7e-e41c-41f2-aeeb-b7ed59250348
begin
	using multiwfs: is_stable, Hcont, Hfilter
	using PlutoUI
	using Plots
	using Base.Meta: parse
end

# ╔═╡ 22b22c60-d257-446d-bc76-305f1e9cdd02
md"""
Writing up the objective properly: we'll say $\text{obj}(V) = \sqrt{\text{obj}_\text{atm}(V) + \text{obj}_\text{NCP}(V) + \text{obj}_\text{noise}(V)}$ where

$$\text{obj}_\text{atm}(V) = \int_{f_\min}^{f_\max} \mathrm{d} f\left|\frac{V}{\phi} (f)\right|^2 \left|\text{PSD atm (f)}\right|$$

$$\text{obj}_\text{NCP}(V) = \int_{f_\min}^{f_\max} \mathrm{d} f\left(\left|\frac{V}{L_\text{fast}} (f)\right|^2 + \left|\frac{V}{L_\text{slow}} (f)\right|^2 \right) \left|\text{PSD NCP (f)}\right|$$

$$\text{obj}_\text{noise}(V) = \int_{f_\min}^{f_\max} \mathrm{d} f\left(\left|\frac{V}{N_\text{fast}} (f)\right|^2 + \left|\frac{V}{N_\text{slow}} (f)\right|^2 \right) \left|\text{PSD noise (f)}\right|$$

and we want to minimize $\text{obj}(X)$ while keeping $\text{obj}(Y)$ below a certain threshold (say 1 rad).
"""

# ╔═╡ 32d82de9-46e7-4a80-94a2-54859ba92ce7
md"""Note: Nslow term should only go up to Nyquist on the slow sensor (50 Hz)"""

# ╔═╡ ab7b9241-99f7-462c-8673-97c82208e5a4


# ╔═╡ e0349f46-24e0-4b40-b578-c98ccddc9916
setting = "optimal" # "optimal" or "lisa" or "slider"

# ╔═╡ c1cc9909-f4dd-47c3-b0d5-4b097135412f
@bind f_cutoff_s Slider(1.0:1.0:25.0)

# ╔═╡ 70a8f1c5-9fd6-4dde-b749-41fad1bb88ac
@bind gain_slow_s Slider(0.01:0.01:2.0)

# ╔═╡ 90676a3e-5808-4d86-9685-c386b6e02ad7
@bind gain_fast_s Slider(0.01:0.01:1.0)

# ╔═╡ 9303f6a9-88a7-4545-bb9f-8656b9b31e74
begin
	f_cutoff, gain_slow, gain_fast = nothing, nothing, nothing
	if setting ==  "optimal"
		f_cutoff, gain_slow, gain_fast = 2.0, 1.74, 0.75
	elseif setting == "lisa"
		f_cutoff, gain_slow, gain_fast = 15.0, 1.4, 0.4
	elseif setting == "slider"
		f_cutoff, gain_slow, gain_fast = f_cutoff_s, gain_slow_s, gain_fast_s
	end
end;

# ╔═╡ 0c8d2a1c-f2e6-4c22-bc0c-c84027bf386b
begin	
	f_loop = 1000.0
	f_noise_crossover = 200.0
	fr = 10 .^ (-2:0.01:log10(f_loop/2))
	sr = 2π .* im .* fr
	sT = sr / f_loop
	ar1_high = ar1_filter(f_cutoff, f_loop, "high")
	no_filter = ZPKFilter(0, 0, 1)
	# f_loop, frame_delay, gain, leak, fpf
	sys_fast = AOSystem(f_loop, 1.0, gain_fast, 0.999, 1, ar1_high)
	sys_slow = AOSystem(f_loop / 10, 0.1, gain_slow, 0.999, 1, no_filter)
	Cfast = sT -> Hcont(sys_fast, sT * f_loop) * Hfilter(sys_fast, sT * f_loop)
	Cslow = sT -> Hcont(sys_slow, sT * f_loop)
	vk_atm = VonKarman(v=10)
	vk_ncp = VonKarman(v=0.1, rms_target=0.8)
	noise_normalization = psd_von_karman(f_noise_crossover, vk_atm)
end;

# ╔═╡ f7c6210d-1e30-489a-85b5-6a892fd6c2f6
	begin
	    plots_contrib = []
	    for (v, sys) in zip(["X", "Y"], [sys_slow, Vector{multiwfs.AOSystem}([sys_slow, sys_fast])])
			tfc = is_stable(sys) ? :green : :red
	        ne = eval(parse("notched_error_$v"))
	        p_v = plot(legend=:bottomleft, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="Closed-loop residual at $v (rad/Hz)", title="\n $v error = $(round(ne(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover), digits=3)) rad \n", titlefontcolor=tfc, ylims=(1e-10, 1e2))
			for errsource in ["atm", "ncp", "noise"]
				err_source_fn = eval(parse("$(errsource)_error_at_f_$v"))
				plot!(fr, err_source_fn.(fr, Cfast, Cslow, 10, Ref(vk_atm), Ref(vk_ncp), noise_normalization, f_loop), label=errsource)
			end
			vline!([f_cutoff], color=:black, ls=:dash, label="HPF cutoff")
	        push!(plots_contrib, p_v)
	    end
	    plot(plots_contrib..., suptitle="$setting: fc = $f_cutoff, gslow = $gain_slow, gfast = $gain_fast")
	end

# ╔═╡ 8fb4f1d2-1ae9-4cbd-9bba-3b02fe13d913
begin
	    plots_etf = []
	    for v in ["X", "Y"]
	        ne = eval(parse("notched_error_$v"))
	        p_v = plot(legend=:bottomright, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", ylims=(1e-10, 1e2), title="$v error = $(round(ne(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover), digits=3)) rad")
	        for (fname, c, s) in zip(["phi_to_$v", "Lfast_to_$v", "Lslow_to_$v", "Nfast_to_$v", "Nslow_to_$v"], [1, 2, 2, 3, 3], [:solid, :solid, :dash, :solid, :dash])
	            f = eval(parse(fname))
	            plot!(fr, abs2.(f.(sT, Cfast, Cslow, 10)), label="|$fname|²", c=c, ls=s)
	        end
	        push!(plots_etf, p_v)
	    end
	    plot(plots_etf...)
	end

# ╔═╡ bad74e0b-f9db-4107-96e0-315d706423f1
begin
	v = Vector{multiwfs.AOSystem}([sys_slow, sys_fast])
	nyquist_plot(v)
end

# ╔═╡ 8131ee09-cb2b-4a0c-8e0b-3780f4899e78
is_stable([sys_slow, sys_fast])

# ╔═╡ 1b292952-30a4-4c18-89b4-c97d29ed4a07
is_stable(Vector{multiwfs.AOSystem}([sys_slow, sys_fast]))

# ╔═╡ Cell order:
# ╠═65e6fa3f-0342-434b-b37d-8a0ec7ddc1f9
# ╠═247e2d7e-e41c-41f2-aeeb-b7ed59250348
# ╟─22b22c60-d257-446d-bc76-305f1e9cdd02
# ╟─32d82de9-46e7-4a80-94a2-54859ba92ce7
# ╠═ab7b9241-99f7-462c-8673-97c82208e5a4
# ╠═0c8d2a1c-f2e6-4c22-bc0c-c84027bf386b
# ╟─9303f6a9-88a7-4545-bb9f-8656b9b31e74
# ╠═e0349f46-24e0-4b40-b578-c98ccddc9916
# ╠═c1cc9909-f4dd-47c3-b0d5-4b097135412f
# ╠═70a8f1c5-9fd6-4dde-b749-41fad1bb88ac
# ╠═90676a3e-5808-4d86-9685-c386b6e02ad7
# ╟─f7c6210d-1e30-489a-85b5-6a892fd6c2f6
# ╟─8fb4f1d2-1ae9-4cbd-9bba-3b02fe13d913
# ╠═bad74e0b-f9db-4107-96e0-315d706423f1
# ╠═8131ee09-cb2b-4a0c-8e0b-3780f4899e78
# ╠═1b292952-30a4-4c18-89b4-c97d29ed4a07
