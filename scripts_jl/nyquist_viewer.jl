### A Pluto.jl notebook ###
# v0.19.43

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

# ╔═╡ b6f41cbe-42eb-11ef-29db-f7b08d9f8ed0
begin
	using Pkg
	Pkg.activate("..")
	using multiwfs
	using Plots
	using PlutoUI
	using Revise
end

# ╔═╡ 0dcefc33-c624-446a-a2f2-7d342bacecf3
@bind f_cutoff Slider(2.0:0.5:100.0)

# ╔═╡ 46c274e6-25f8-483a-9a06-ccb7d0fb3b2e
@bind delay Slider(0.0:0.05:1.0)

# ╔═╡ 30ed7005-7e8c-465d-a0c8-4219ad7239a7
delay

# ╔═╡ 348b8b94-e344-48bf-8e66-09947d1e91ca
sys = AOSystem(1000.0, delay, 0.01, 0.999, 10, "high", f_cutoff);

# ╔═╡ a6f1b032-4dba-4cd9-ada9-7be9e8829bff
round(search_gain!(sys), digits=2)

# ╔═╡ 264e2435-10de-4217-8bf7-0d65fd837723
nyquist_plot(sys)

# ╔═╡ Cell order:
# ╠═b6f41cbe-42eb-11ef-29db-f7b08d9f8ed0
# ╠═0dcefc33-c624-446a-a2f2-7d342bacecf3
# ╠═46c274e6-25f8-483a-9a06-ccb7d0fb3b2e
# ╠═30ed7005-7e8c-465d-a0c8-4219ad7239a7
# ╠═348b8b94-e344-48bf-8e66-09947d1e91ca
# ╠═a6f1b032-4dba-4cd9-ada9-7be9e8829bff
# ╠═264e2435-10de-4217-8bf7-0d65fd837723
