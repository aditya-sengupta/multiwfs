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

# ╔═╡ 62aca278-ff69-11ef-1e27-13496cd4b244
begin
	using Pkg
	Pkg.activate("..")
	using Plots
	using multiwfs
	using PlutoUI
end

# ╔═╡ 78218df1-3b7f-4908-91ea-1d25842012e7
begin
	r0 = 0.1031
	r0_ncp = 0.6
	crossover_freq = 500.0
	vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
	vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
	R = 10
	f_loop = 1000.0
end

# ╔═╡ 8ba46896-2905-4818-a738-bb0d7edddf3b
begin
	A_ar1 = [0.999 0; 1 0]
    L = A_DM(2)
    Ã = block_diag(L, A_ar1)
    C̃ = [0 -1 0 1]
    D̃ = [1 0 0 0]' 
    B = [0; 0; 1; 0]
    Pw = hcat(1...)
    W = B * Pw * B'
    V = hcat(0.2...)
    K̃ = kalman_gain(Ã, C̃, W, V)
    Vv = [0 -1 0 1]
    Q = Vv' * Vv
    Rlqg = zeros(1,1)
    L = lqr_gain(Ã, D̃, Q, Rlqg)
    lqgic_controller = LQG(Ã, D̃, C̃, K̃, L, 1/f_loop)
end

# ╔═╡ 8373c525-b700-42ee-ae8d-c699912d0dff
sT_cr = 5.0 * 2π * im / 1000.0

# ╔═╡ 2c8883bf-0ad8-46c4-897c-13df6bc1cfc6
@bind gain_slow Slider(0.0:0.01:2.0)

# ╔═╡ 1d9ec936-594b-404a-a6b2-a83acf468298
@bind gain_fast Slider(0.01:0.01:1.0)

# ╔═╡ b74c1d52-7b9b-43e5-aa5c-cae6a9d99ff0
@bind cutoff_freq Slider(0.0:0.5:100.0)

# ╔═╡ 8017016c-5a9f-4dca-a515-3e5f36a8788c
gain_slow, gain_fast, cutoff_freq

# ╔═╡ c537df79-20f1-4407-82d5-6528d3794aea
begin
	fmc = 1.0
	no_filter = ZPKFilter(0, 0, 1)
	slow_controller = FilteredIntegrator(gain_slow, 0.999, ZPKFilter(0, 0, 1), R/f_loop)
	fast_controller_unfiltered = FilteredIntegrator(gain_fast, 0.999, no_filter, 1/f_loop)
	fast_controller_filtered = FilteredIntegrator(gain_fast, 0.999, ar1_filter(cutoff_freq, f_loop, "high"), 1/f_loop)
	sim_filtered = Simulation(f_loop, fast_controller_filtered, slow_controller, R, vk_atm, vk_ncp, crossover_freq, f_min_cost=fmc)
	sim_unfiltered = Simulation(f_loop, fast_controller_unfiltered, slow_controller, R, vk_atm, vk_ncp, crossover_freq, f_min_cost=fmc)
	sim_lqgic = Simulation(f_loop, lqgic_controller, slow_controller, R, vk_atm, vk_ncp, crossover_freq, f_min_cost=fmc)
end

# ╔═╡ 6d288c29-2d37-4d6d-9e49-c5f25fe8ccf9
abs2.(phi_to_X(sT_cr, sim_filtered))

# ╔═╡ 12582d67-1908-4910-85a3-cfe53658e150
abs2(Lfast_to_X(sT_cr, sim_filtered)) + abs2(Lslow_to_X(sT_cr, sim_filtered))

# ╔═╡ 00ee7711-6034-4dd6-95e7-592c13957da7
abs2(Nfast_to_X(sT_cr, sim_filtered)) + abs2(Nslow_to_X(sT_cr, sim_filtered))

# ╔═╡ e21066e6-9862-4c5e-8786-083813c9bbd9
ne_filt = notched_error_X(sim_filtered)

# ╔═╡ f2db9b1e-6664-452b-bb16-e75b4c1bc9e0
ne_unfilt = notched_error_X(sim_unfiltered)

# ╔═╡ 498db802-0d10-4620-94c8-e0cb6233306b
function ne_from_params(gain_slow, gain_fast, cutoff_freq)
	slow_controller = FilteredIntegrator(gain_slow, 0.999, no_filter, R/f_loop)
	fast_controller_filtered = FilteredIntegrator(gain_fast, 0.999, ar1_filter(cutoff_freq, f_loop, "high"), 1/f_loop)
	sim_filtered = Simulation(f_loop, fast_controller_filtered, slow_controller, R, vk_atm, vk_ncp, crossover_freq, f_min_cost=fmc)
	return notched_error_X(sim_filtered)
end

# ╔═╡ f66ab79d-926f-4c99-9acc-e3e056e3e655
lisa_notched_error = ne_from_params(1.4, 0.4, 15.0)

# ╔═╡ 0ca5084c-2982-4ea1-b18d-f2cbfd472126
ne_from_params(gain_slow, gain_fast, cutoff_freq)

# ╔═╡ 807ae7e1-ff1a-4b0f-b093-1a8e1fcac8db
(ne_from_params(gain_slow, gain_fast, cutoff_freq) - lisa_notched_error) / lisa_notched_error * 100

# ╔═╡ c26b76ba-7b45-4731-8382-b329719dff14
ne_from_params(0.0, 0.4, 15.0)

# ╔═╡ 1cbe4648-8c48-47be-bbf8-79a8ab2308bc
p_filtered = begin
	plot(sim_filtered.fr, atm_error_at_f_X.(sim_filtered.fr, Ref(sim_filtered)), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", label="atm", title="Filtered, X err = $(round(ne_filt, digits=3))", ylim=(1e-6, 1e2), xticks=exp10.(-3:2))
	plot!(sim_filtered.fr, ncp_error_at_f_X.(sim_filtered.fr, Ref(sim_filtered)), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", label="ncp")
	plot!(sim_filtered.fr, noise_error_at_f_X.(sim_filtered.fr, Ref(sim_filtered)), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", label="noise")
	vline!([fmc], ls=:dash, color=:black, label="Cost cutoff")
	vline!([cutoff_freq], ls=:dash, color=:purple, label="HPF cutoff")
	vline!([zero_db_bandwidth(sim_filtered)], ls=:dash, color=:teal, label="Bandwidth")
end;

# ╔═╡ 494e929c-9335-4b0e-acd5-d9df5650dff1
p_filtered

# ╔═╡ d468cc84-d776-4a05-9851-77a95deae471
p_unfiltered = begin
	plot(sim_unfiltered.fr, atm_error_at_f_X.(sim_unfiltered.fr, Ref(sim_unfiltered)), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", label="atm", title="Unfiltered, X err = $(round(ne_unfilt, digits=3))", ylim=(1e-6, 1e2), xticks=exp10.(-3:2))
	plot!(sim_unfiltered.fr, ncp_error_at_f_X.(sim_unfiltered.fr, Ref(sim_unfiltered)), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", label="ncp")
	plot!(sim_unfiltered.fr, noise_error_at_f_X.(sim_unfiltered.fr, Ref(sim_unfiltered)), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", label="noise")
	vline!([fmc], ls=:dash, color=:black, label="Cost cutoff")
end;

# ╔═╡ 1eb92812-2945-4b68-9c54-ecead942e9a7
nyquist_plot(sim_filtered)

# ╔═╡ e0a0948c-9d5c-4703-8280-fac7370d19d9


# ╔═╡ Cell order:
# ╠═62aca278-ff69-11ef-1e27-13496cd4b244
# ╠═78218df1-3b7f-4908-91ea-1d25842012e7
# ╠═8ba46896-2905-4818-a738-bb0d7edddf3b
# ╠═8373c525-b700-42ee-ae8d-c699912d0dff
# ╠═6d288c29-2d37-4d6d-9e49-c5f25fe8ccf9
# ╠═12582d67-1908-4910-85a3-cfe53658e150
# ╠═00ee7711-6034-4dd6-95e7-592c13957da7
# ╠═f66ab79d-926f-4c99-9acc-e3e056e3e655
# ╠═0ca5084c-2982-4ea1-b18d-f2cbfd472126
# ╠═807ae7e1-ff1a-4b0f-b093-1a8e1fcac8db
# ╠═2c8883bf-0ad8-46c4-897c-13df6bc1cfc6
# ╠═1d9ec936-594b-404a-a6b2-a83acf468298
# ╠═b74c1d52-7b9b-43e5-aa5c-cae6a9d99ff0
# ╠═8017016c-5a9f-4dca-a515-3e5f36a8788c
# ╠═c26b76ba-7b45-4731-8382-b329719dff14
# ╠═494e929c-9335-4b0e-acd5-d9df5650dff1
# ╠═c537df79-20f1-4407-82d5-6528d3794aea
# ╠═e21066e6-9862-4c5e-8786-083813c9bbd9
# ╠═f2db9b1e-6664-452b-bb16-e75b4c1bc9e0
# ╠═498db802-0d10-4620-94c8-e0cb6233306b
# ╠═1cbe4648-8c48-47be-bbf8-79a8ab2308bc
# ╠═d468cc84-d776-4a05-9851-77a95deae471
# ╠═1eb92812-2945-4b68-9c54-ecead942e9a7
# ╠═e0a0948c-9d5c-4703-8280-fac7370d19d9
