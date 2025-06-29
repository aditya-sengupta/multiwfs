using multiwfs
using NPZ
using Plots
using QuadGK
using JLD2
using ProgressMeter
using Statistics
using multiwfs: fast_ncp_error_at_f_X, slow_ncp_error_at_f_X, fast_noise_error_at_f_X, slow_noise_error_at_f_X

db_Nfreqpoints = 1000
Nfreqpoints = db_Nfreqpoints ÷ 2
k = Nfreqpoints + 1
freq = range(0.0, 500.0, k)[2:end]
f_crossover = 200.0

begin
    optpars = load("data/all_opt.jld2")
    xerrs = Dict()
    r0_ncps = vcat(0.6:0.2:2.0, 3:6)    
    for r0_ncp in r0_ncps
        vk_ncp = VonKarman(0.001, 0.25 * (r0_ncp)^(-5/3))
        sim_one = simgen_ichpf(optpars["($r0_ncp, $f_crossover)"][1][1], 0.0, 0.0, vk_ncp, f_crossover; leak=0.999)
        sim_two = simgen_ichpf(optpars["($r0_ncp, $f_crossover)"][3]..., vk_ncp, f_crossover; leak=0.999)
        xerrs["($r0_ncp, $f_crossover) one"] = notched_error_X(sim_one)
        xerrs["($r0_ncp, $f_crossover) two"] = notched_error_X(sim_two)
    end
    
    p = Plots.plot(title="fcross = $f_crossover Hz", xlabel="r₀ NCP (m)", ylabel="CL error at X (rad)", xscale=:log10, xticks=([0.6, 1.0, 1.4, 2.0, 3.0, 4.0, 6.0], [0.6, 1.0, 1.4, 2.0, 3.0, 4.0, 6.0]))
    Plots.plot!(r0_ncps, [xerrs["($r0_ncp, $f_crossover) one"] for r0_ncp in r0_ncps], color=1, msw=0, legend=nothing)
    Plots.plot!(r0_ncps, [xerrs["($r0_ncp, $f_crossover) two"] for r0_ncp in r0_ncps], color=2, msw=0, legend=nothing)

    Plots.plot(p)
end

