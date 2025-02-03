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

# ╔═╡ 94bf96fd-0faa-4b11-b1a3-0e27b81a4258
begin
	using Pkg
	Pkg.activate("..")
	using multiwfs
	using multiwfs: Hcont, Hfilter
	using PlutoUI
	using Base.Meta: parse
	using Plots
end

# ╔═╡ 22b22c60-d257-446d-bc76-305f1e9cdd02
md"""
Writing up the objective properly: we'll say $\text{obj}(V) = \sqrt{\text{obj}_\text{atm}(V) + \text{obj}_\text{NCP}(V) + \text{obj}_\text{noise}(V)}$ where

$$\text{obj}_\text{atm}(V) = \int \mathrm{d} f\left|\frac{V}{\phi}\right|^2 \left|\text{PSD atm}\right|$$

$$\text{obj}_\text{NCP}(V) = \int \mathrm{d} f\left(\left|\frac{V}{L_\text{fast}}\right|^2 + \left|\frac{V}{L_\text{slow}}\right|^2 \right) \left|\text{PSD NCP}\right|$$

$$\text{obj}_\text{noise}(V) = \int \mathrm{d} f\left(\left|\frac{V}{N_\text{fast}}\right|^2 + \left|\frac{V}{N_\text{slow}}\right|^2 \right) \left|\text{PSD noise}\right|$$

and we want to minimize $\text{obj}(X)$ while keeping $\text{obj}(Y)$ below a certain threshold (say 1 rad).
"""

# ╔═╡ c1cc9909-f4dd-47c3-b0d5-4b097135412f
@bind f_cutoff Slider(1.0:1.0:100.0)

# ╔═╡ 70a8f1c5-9fd6-4dde-b749-41fad1bb88ac
@bind gain_slow Slider(0.01:0.01:2.0)

# ╔═╡ 90676a3e-5808-4d86-9685-c386b6e02ad7
@bind gain_fast Slider(0.01:0.01:1.0)

# ╔═╡ a920f762-4381-4002-8535-ee08885021dc
begin
	f_loop = 1000.0
	f_noise_crossover = 200.0
	fr = 10 .^ (-2:0.01:log10(f_loop/2))
	sr = 2π .* im .* fr
	sT = sr / f_loop
	ar1_high = ar1_filter(f_cutoff, f_loop, "high")
	no_filter = ZPKFilter(0, 0, 1)
	# f_loop, frame_delay, gain, leak, fpf
	sys_high = AOSystem(f_loop, 1.0, gain_slow, 0.999, 1, ar1_high)
	sys_low = AOSystem(f_loop, 1.0, gain_fast, 0.999, 1, no_filter)
	Cfast = sT -> Hcont(sys_high, sT * f_loop) * Hfilter(sys_high, sT * f_loop)
	Cslow = sT -> Hcont(sys_low, sT * f_loop)
	vk_atm = VonKarman(v=10)
	vk_ncp = VonKarman(v=0.1, rms_target=0.8)
	noise_normalization = psd_von_karman(200.0, vk_atm)
end;

# ╔═╡ 1062b63d-13b5-4ede-943f-02d0d958e9f2
begin
    plots = []
    for v in ["X", "Y"]
        ne = eval(parse("notched_error_$v"))
        p_v = plot(legend=:bottomleft, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="Closed-loop residual at $v (rad)", title="$v error = $(round(ne(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover), digits=3)) rad")
		for errsource in ["atm", "ncp", "noise"]
			err_source_fn = eval(parse("$(errsource)_error_at_f_$v"))
			plot!(fr, err_source_fn.(fr, Cfast, Cslow, 10, Ref(vk_atm), Ref(vk_ncp), noise_normalization, f_loop), label=errsource)
		end
		vline!([f_cutoff], color=:black, ls=:dash, label="HPF cutoff")
        push!(plots, p_v)
    end
    plot(plots...)
end

# ╔═╡ Cell order:
# ╟─94bf96fd-0faa-4b11-b1a3-0e27b81a4258
# ╟─22b22c60-d257-446d-bc76-305f1e9cdd02
# ╠═c1cc9909-f4dd-47c3-b0d5-4b097135412f
# ╠═70a8f1c5-9fd6-4dde-b749-41fad1bb88ac
# ╠═90676a3e-5808-4d86-9685-c386b6e02ad7
# ╠═1062b63d-13b5-4ede-943f-02d0d958e9f2
# ╟─a920f762-4381-4002-8535-ee08885021dc
