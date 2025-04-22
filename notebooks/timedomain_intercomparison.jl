### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# ╔═╡ 931ee824-1e5e-11f0-3975-b32e2429c3cc
begin
	using Pkg
	Pkg.activate("..")
	using NPZ
	using PlutoLinks
	PlutoLinks.@revise using multiwfs
end

# ╔═╡ 4b64ae58-917d-4d1f-ab6f-d9c31cd07a5c
(atm=npzread("../data/olt/olt_atm.npy"))

# ╔═╡ Cell order:
# ╠═931ee824-1e5e-11f0-3975-b32e2429c3cc
# ╠═4b64ae58-917d-4d1f-ab6f-d9c31cd07a5c
