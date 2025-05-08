using multiwfs
using JLD2
using DataFrames
using CSV
using DataStructures: OrderedDict

optpars_table = OrderedDict(
    "r0_ncp" => [],
    "fcross" => [],
    "gain_onewfs" => [],
    "error_onewfs" => [],
    "gain_fast_twowfs" => [],
    "gain_slow_twowfs" => [],
    "f_cutoff_twowfs" => [],
    "alpha_twowfs" => [],
    "error_twowfs" => [],
    "gain_fast_twowfs_unlim" => [],
    "gain_slow_twowfs_unlim" => [],
    "f_cutoff_twowfs_unlim" => [],
    "alpha_twowfs_unlim" => [],
    "error_twowfs_unlim" => [],
)

optpars_ic = load("data/ic_opt_for_paper.jld2")
optpars_ichpf = load("data/ic_opt_lisa_20250430.jld2")
for r0_ncp in [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, 6.0]
    vk_ncp = VonKarman(0.001, 0.25 * (r0_ncp)^(-5/3))
    for fcross in [50.0, 100.0, 200.0, 500.0]
        pars_onewfs = optpars_ic["($r0_ncp, $fcross)"][1]
        push!(optpars_table["r0_ncp"], r0_ncp)
        push!(optpars_table["fcross"], fcross)
        push!(optpars_table["gain_onewfs"], pars_onewfs[1])
        sim_onewfs = simgen_ichpf(pars_onewfs[1], 0.0, 0.0, vk_ncp, fcross)
        push!(optpars_table["error_onewfs"], notched_error_X(sim_onewfs))

        pars_twowfs = optpars_ichpf["($r0_ncp, $fcross)"][2]
        pars_twowfs_unlim = optpars_ichpf["($r0_ncp, $fcross)"][1]
        for (label, partwo) in zip(["", "_unlim"], [pars_twowfs, pars_twowfs_unlim])
            sim = simgen_ichpf(partwo..., vk_ncp, fcross)
            push!(optpars_table["gain_fast_twowfs" * label], partwo[1])
            push!(optpars_table["gain_slow_twowfs" * label], partwo[2])
            push!(optpars_table["f_cutoff_twowfs" * label], partwo[3])
            α = exp(-2π * partwo[3] / sim.f_loop)
            push!(optpars_table["alpha_twowfs" * label], α)
            push!(optpars_table["error_twowfs" * label], notched_error_X(sim))
        end
    end
end

CSV.write("data/optpars_20250501.csv", DataFrame(optpars_table))