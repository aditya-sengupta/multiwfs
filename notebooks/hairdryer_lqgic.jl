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
	using PlutoLinks
	using PlutoUI
end

# ╔═╡ 32cf0c90-9f03-48e3-81ea-5aea6ccef557
PlutoLinks.@revise using multiwfs

# ╔═╡ 78218df1-3b7f-4908-91ea-1d25842012e7
begin
	r0 = 0.1031
	r0_ncp = 0.6
	crossover_freq = 500.0
	vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
	vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
	R = 10
	f_loop = 1000.0
	fmc = 1.0
	no_filter = ZPKFilter(0, 0, 1)
end

# ╔═╡ d5efb157-c3b4-4caa-a77c-715f323e8fa6
@bind gain_slow Slider(0.0:0.01:2.0)

# ╔═╡ d787d852-8dd6-4b7e-9817-cbd30164e32c
@bind log_lqg_noise Slider(-4:0.01:2)

# ╔═╡ e5b0e815-7fb0-47a7-b2e2-8f9ea47db376
@bind num_nines_slow Slider(1:0.01:5)

# ╔═╡ 9fcb915b-82b4-45f2-b011-41396b6f2b6e
@bind num_nines_fast Slider(1:0.01:5)

# ╔═╡ 7d1b40b7-87ea-498e-8281-d0f107b80776
leak_slow, leak_fast = 1 - exp10(-num_nines_slow), 1 - exp10(-num_nines_fast);

# ╔═╡ df6693b9-357d-42bb-a67a-bde67d31e8ec
@bind f_cutoff Slider(0.0:0.1:100.0)

# ╔═╡ 8ba46896-2905-4818-a738-bb0d7edddf3b
begin
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
    lqgic_controller = FilteredLQG(LQG(Ã, D̃, C̃, K̃, L, 1/f_loop), ar1_filter(f_cutoff, f_loop, "high"))
	lisa_fast_controller = FilteredIntegrator(0.4, 0.999, ar1_filter(15.0, 1000.0, "high"), 1e-3)
	lisa_slow_controller = FilteredIntegrator(1.4, 0.999, no_filter, 1e-2)
end

# ╔═╡ 2c7aa02a-013a-42a1-bc2f-102cf746bba0
sim_lisa = Simulation(f_loop, lisa_fast_controller, lisa_slow_controller, R, vk_atm, vk_ncp, crossover_freq, f_min_cost=fmc)

# ╔═╡ d558fd99-102a-488a-afc6-9c924a7009b2
lisa_notched_error = notched_error_X(sim_lisa);

# ╔═╡ 3b006068-29cb-4da2-a182-05a27030f681
gain_slow, log_lqg_noise, num_nines_slow, num_nines_fast, f_cutoff

# ╔═╡ c537df79-20f1-4407-82d5-6528d3794aea
begin
	slow_controller = FilteredIntegrator(gain_slow, leak_slow, no_filter, R/f_loop)
	sim_lqgic = Simulation(f_loop, lqgic_controller, slow_controller, R, vk_atm, vk_ncp, crossover_freq, f_min_cost=fmc)
end;

# ╔═╡ 46836717-743b-4a6e-8992-a10ddf7be688
(lisa_notched_error - notched_error_X(sim_lqgic)) / lisa_notched_error

# ╔═╡ 8dba1925-7692-4675-8229-91f60c343e15
plot(
	plot_integrands(sim_lqgic; legend=:topright, xticks=exp10.(-3:2)),
	nyquist_plot(sim_lqgic, legend=nothing, title="X err = $(round(notched_error_X(sim_lqgic), digits=3)) rad")
)

# ╔═╡ 24c99b31-e956-4dc5-9bde-061469aae29d
plot(
	plot_integrands(sim_lisa; legend=:topright, xticks=exp10.(-3:2)),
	nyquist_plot(sim_lisa, legend=nothing, title="X err = $(round(notched_error_X(sim_lisa), digits=3)) rad")
)

# ╔═╡ Cell order:
# ╠═62aca278-ff69-11ef-1e27-13496cd4b244
# ╠═32cf0c90-9f03-48e3-81ea-5aea6ccef557
# ╠═78218df1-3b7f-4908-91ea-1d25842012e7
# ╠═8ba46896-2905-4818-a738-bb0d7edddf3b
# ╠═2c7aa02a-013a-42a1-bc2f-102cf746bba0
# ╠═d558fd99-102a-488a-afc6-9c924a7009b2
# ╟─7d1b40b7-87ea-498e-8281-d0f107b80776
# ╠═3b006068-29cb-4da2-a182-05a27030f681
# ╠═d5efb157-c3b4-4caa-a77c-715f323e8fa6
# ╠═d787d852-8dd6-4b7e-9817-cbd30164e32c
# ╠═e5b0e815-7fb0-47a7-b2e2-8f9ea47db376
# ╠═9fcb915b-82b4-45f2-b011-41396b6f2b6e
# ╠═df6693b9-357d-42bb-a67a-bde67d31e8ec
# ╠═c537df79-20f1-4407-82d5-6528d3794aea
# ╠═46836717-743b-4a6e-8992-a10ddf7be688
# ╟─8dba1925-7692-4675-8229-91f60c343e15
# ╠═24c99b31-e956-4dc5-9bde-061469aae29d
