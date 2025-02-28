### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 62bca1df-4593-493b-99b9-3a0382f18bc6
begin
	using Pkg
	Pkg.activate("..")
	using PlutoLinks: @revise	
end

# ╔═╡ a057dd55-5181-46e5-b860-b6ccb2974179
begin
	@revise using multiwfs
	using multiwfs: Hcont, Hfilter, is_stable
	using Plots
	using Base: product
	using Base.Meta: parse
	using Base.Threads: @threads
	using ProgressMeter: @showprogress
end

# ╔═╡ c282dcbe-c6e0-4c8f-b2b9-f412a0231702
begin
	f_loop = 1000.0
	fr = 10 .^ (-2:0.01:log10(f_loop/2))
	sr = 2π .* im .* fr
	sT = sr / f_loop
	no_filter = ZPKFilter(0, 0, 1)
	vk_atm = VonKarman(v=10)
	vk_ncp = VonKarman(v=0.01, rms_target=8.0)
	f_noise_crossover = 200.0
	noise_normalization = psd_von_karman(200.0, vk_atm)
	sys_slow_orig = AOSystem(f_loop / 10, 0.1, 1.77, 0.999, 1, no_filter)
end;

# ╔═╡ fff0d9d8-e73b-487c-8f1a-520d6526d7ef
sys_fast = AOSystem(f_loop, 1.0, 0.4, 0.999, 1, ar1_filter(1.0, f_loop, "high"))

# ╔═╡ da49cc8f-10f6-4527-95b6-86276d628592
nyquist_plot(sys_fast)

# ╔═╡ 0fbfdc7e-79a3-49bf-851a-3c5f8be9bb17
md"""
Let's make the unfiltered slow system in LQG form and make sure we get the same open-loop TF
"""

# ╔═╡ e75fdaa1-068a-40be-a0f2-55ab8e9441de
begin
	α = 0.995
	β = -0.01
	Vc = 1
	A_ar1 = [α β; 1 0]
    L = A_DM(2)
    Ã = block_diag(L, A_ar1)
    C̃ = [0 -1 0 1]
    D̃ = [1 0 0 0]' 
    B = [0; 0; 1; 0]
    Pw = hcat(1...)
    W = B * Pw * B'
    V = hcat(Vc...)
    K̃ = kalman_gain(Ã, C̃, W, V)
    Vv = [0 -1 0 1]
    Q = Vv' * Vv
    R = zeros(1,1)
    L = lqr_gain(Ã, D̃, Q, R)
    lqg = LQG(Ã, D̃, C̃, K̃, L)
	sys_slow = AOSystem(f_loop / 10, 0.1, 1.77, 0.999, 1, lqg)
end;

# ╔═╡ 2c0446ff-244c-4a2c-beb7-3c58e09191e2
begin
	nyquist_plot(Vector{AOSystem}([sys_fast, sys_slow_orig]))
	nyquist_plot!(Vector{AOSystem}([sys_fast, sys_slow]))
end

# ╔═╡ 306a445e-e4ed-48f3-bfee-4baf60208b4e
Cfast = sT -> Hcont(sys_fast, sT * f_loop) * Hfilter(sys_fast, sT * f_loop)

# ╔═╡ fda4f256-bd01-403a-becc-bc3cc6b154e0
Cslow = sT -> Hcont(sys_slow, sT * f_loop)

# ╔═╡ 4ee50ab9-d1ea-43b7-a88f-4845314ede7e
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

# ╔═╡ 6a8a649d-5f33-4c17-8f22-a24a251c6403
begin
	 plots_contrib = []
	stable = is_stable(Vector{multiwfs.AOSystem}([sys_slow, sys_fast]))
	tfc = stable ? :green : :red
	errX, errY = notched_error_X(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover), notched_error_Y(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover)
	p_v = plot(legend=:bottomleft, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="Closed-loop residual (rad/Hz)", title="X error = $(round(errX, digits=3)) rad, Y error = $(round(errY, digits=3)) rad", titlefontcolor=tfc, ylims=(1e-10, 1e2))
	for errsource in ["atm_error_at_f_X", "ncp_error_at_f_X", "ncp_error_at_f_Y", "noise_error_at_f_X"]
		err_source_fn = eval(parse(errsource))
		plot!(fr, err_source_fn.(fr, Cfast, Cslow, 10, Ref(vk_atm), Ref(vk_ncp), noise_normalization, f_loop), label=errsource)
	end
	p_v
end

# ╔═╡ 62aaaaa2-d50f-44bd-8e58-e85c321046fe
begin
	plot(fr, psd_von_karman.(fr, Ref(vk_atm)), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)", label="Atmosphere input PSD")
	plot!(fr, psd_von_karman.(fr, Ref(vk_ncp)), label="NCP input PSD")
	hline!([psd_von_karman(200, vk_atm)], label="Noise", legend=:topright)
end

# ╔═╡ 16780668-02e8-46b3-9a78-37e8d367258a
vk_atm.f₀

# ╔═╡ 6a1d8dbb-2490-4cd2-84cd-032bae515200
vk_ncp.f₀

# ╔═╡ Cell order:
# ╠═62bca1df-4593-493b-99b9-3a0382f18bc6
# ╠═a057dd55-5181-46e5-b860-b6ccb2974179
# ╠═c282dcbe-c6e0-4c8f-b2b9-f412a0231702
# ╠═fff0d9d8-e73b-487c-8f1a-520d6526d7ef
# ╠═da49cc8f-10f6-4527-95b6-86276d628592
# ╟─0fbfdc7e-79a3-49bf-851a-3c5f8be9bb17
# ╟─4ee50ab9-d1ea-43b7-a88f-4845314ede7e
# ╠═e75fdaa1-068a-40be-a0f2-55ab8e9441de
# ╟─6a8a649d-5f33-4c17-8f22-a24a251c6403
# ╠═2c0446ff-244c-4a2c-beb7-3c58e09191e2
# ╠═306a445e-e4ed-48f3-bfee-4baf60208b4e
# ╠═fda4f256-bd01-403a-becc-bc3cc6b154e0
# ╠═62aaaaa2-d50f-44bd-8e58-e85c321046fe
# ╠═16780668-02e8-46b3-9a78-37e8d367258a
# ╠═6a1d8dbb-2490-4cd2-84cd-032bae515200
