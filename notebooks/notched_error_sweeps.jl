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

# ╔═╡ 630a26e6-78c4-46a5-a560-124fe869e60b
begin
	using Pkg
	Pkg.activate("..")
	using multiwfs
	using PlutoUI
end

# ╔═╡ 41f1a889-07b9-475b-919e-a21aa9fc60ca
using multiwfs: Hcont, Hfilter

# ╔═╡ 9cf183c3-8679-4625-a55a-a60eb0f44bc5
using Plots

# ╔═╡ 695c553c-d32f-4f51-a787-6f41d969e53e
using Base.Meta: parse

# min |X / phi|^2 |PSD atm| + (|X / Lfast|^2 + |X / Lslow|^2) |PSD NCP| 
# subject to |Y / phi|^2 |PSD atm| + (|Y / Lfast|^2 + |Y / Lslow|^2) |PSD NCP| < threshold

# ╔═╡ d1f3ab96-4caf-4e76-aff0-3587e23ac9dd
f_loop = 1000.0

# ╔═╡ 06059ebc-c1ab-4e33-801b-0114ab76c9b1
fr = 10 .^ (-2:0.01:log10(f_loop/2))

# ╔═╡ 02b78cb4-77a4-43ef-91b3-8bb1a364d630
sr = 2π .* im .* fr / f_loop

# ╔═╡ ba52c496-e9a3-4c5b-8214-6f714cd80a4d
zr = exp.(sr)

# ╔═╡ b129971a-75b0-4435-beae-2e2dc46d947b
no_filter = ZPKFilter(0, 0, 1)
# f_loop, frame_delay, gain, leak, fpf

# ╔═╡ 572461a5-d939-4073-9aac-f6a6c1610429
vk_atm = VonKarman(v=10)

# ╔═╡ ac0b78e6-bb39-4710-8f1d-19c2473b0b6d
vk_ncp = VonKarman(v=0.1, rms_target=0.8)

# first thing to look at, do my ETFs here match up with Lisa's?

# ╔═╡ b37d9fe2-840f-4f3c-bbbf-2258440c5cce
@bind gain_fast Slider(0.0:0.01:1.0)

# ╔═╡ ad911c62-a4c9-4da6-8565-81e97b8f1eb6
@bind gain_slow Slider(0.0:0.01:1.5)

# ╔═╡ 80a52ee7-de67-4563-b9cd-8db7d9ab7317
sys_low = AOSystem(f_loop / 10, 0.0, gain_slow, 0.999, 1, no_filter)

# ╔═╡ 44a4cc63-8651-47ec-b988-9d8ae3610bf4
Cslow = z -> Hcont(sys_low, log(z))

# ╔═╡ 0491dc7b-55ac-4e10-89e4-0074b30035c6
@bind hpf_cutoff Slider(0.0:1.0:100.0)

# ╔═╡ dd0550b0-92d6-47df-bc23-3bea411c8fd0
ar1_high = ar1_filter(hpf_cutoff, f_loop, "high")

# ╔═╡ cee8e06a-5601-423c-bf28-ea4bb861882b
sys_high = AOSystem(f_loop, 1.0, gain_fast, 0.999, 10, ar1_high)

# ╔═╡ 7b64c612-ff7c-4c34-a0a8-cfb663813130
Cfast = z -> Hcont(sys_high, log(z)) * Hfilter(sys_high, log(z))

# ╔═╡ 03c23a8f-8277-4c14-b42c-6aad02de97d1
notched_error_X(Cfast, Cslow, 10, vk_atm, vk_ncp)

# ╔═╡ 57ac0000-c176-4527-a140-d4ce1b726cb0
notched_error_Y(Cfast, Cslow, 10, vk_atm, vk_ncp)

# ╔═╡ a745114d-118e-4c7d-89a1-20f271dd63f8
begin
    p_x = plot(legend=:left, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", ylims=(1e-10, 1e1), title="X error = $(round(notched_error_X(Cfast, Cslow, 10, vk_atm, vk_ncp), digits=3)) rad")
    for fname in ["phi_to_X", "Lfast_to_X", "Lslow_to_X"]
        f = eval(parse(fname))
        plot!(fr, abs2.(f.(zr, Cfast, Cslow, 10)), label="|$fname|²",)
    end
    p_y = plot(legend=:left, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", ylims=(1e-10, 1e1), title="Y error = $(round(notched_error_Y(Cfast, Cslow, 10, vk_atm, vk_ncp), digits=3)) rad")
    for fname in ["phi_to_Y", "Lfast_to_Y", "Lslow_to_Y"]
        f = eval(parse(fname))
        plot!(fr, abs2.(f.(zr, Cfast, Cslow, 10)), label="|$fname|²",)
    end
    plot(p_x, p_y)
end

# ╔═╡ Cell order:
# ╠═630a26e6-78c4-46a5-a560-124fe869e60b
# ╠═41f1a889-07b9-475b-919e-a21aa9fc60ca
# ╠═9cf183c3-8679-4625-a55a-a60eb0f44bc5
# ╠═695c553c-d32f-4f51-a787-6f41d969e53e
# ╠═d1f3ab96-4caf-4e76-aff0-3587e23ac9dd
# ╠═06059ebc-c1ab-4e33-801b-0114ab76c9b1
# ╠═02b78cb4-77a4-43ef-91b3-8bb1a364d630
# ╠═ba52c496-e9a3-4c5b-8214-6f714cd80a4d
# ╠═dd0550b0-92d6-47df-bc23-3bea411c8fd0
# ╠═b129971a-75b0-4435-beae-2e2dc46d947b
# ╠═cee8e06a-5601-423c-bf28-ea4bb861882b
# ╠═80a52ee7-de67-4563-b9cd-8db7d9ab7317
# ╠═7b64c612-ff7c-4c34-a0a8-cfb663813130
# ╠═44a4cc63-8651-47ec-b988-9d8ae3610bf4
# ╠═572461a5-d939-4073-9aac-f6a6c1610429
# ╠═ac0b78e6-bb39-4710-8f1d-19c2473b0b6d
# ╠═03c23a8f-8277-4c14-b42c-6aad02de97d1
# ╠═57ac0000-c176-4527-a140-d4ce1b726cb0
# ╠═b37d9fe2-840f-4f3c-bbbf-2258440c5cce
# ╠═ad911c62-a4c9-4da6-8565-81e97b8f1eb6
# ╠═0491dc7b-55ac-4e10-89e4-0074b30035c6
# ╠═a745114d-118e-4c7d-89a1-20f271dd63f8
