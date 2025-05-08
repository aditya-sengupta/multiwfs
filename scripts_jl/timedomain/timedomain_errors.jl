using multiwfs
using NPZ
using Plots
using QuadGK
using JLD2
using ProgressMeter
using Statistics
using multiwfs: fast_ncp_error_at_f_X, slow_ncp_error_at_f_X, fast_noise_error_at_f_X, slow_noise_error_at_f_X

sim = simgen_ichpf(0.4, 1.4, 15.0, VonKarman(0.001, 0.25), 50.0; leak=0.999)

db_Nfreqpoints = 1000
Nfreqpoints = db_Nfreqpoints ÷ 2
k = Nfreqpoints + 1
freq = range(0.0, 500.0, k)[2:end]

begin
    p = Plots.plot()
    for fcross in [50.0, 200.0]
        sim = simgen_ichpf(0.4, 1.4, 15.0, VonKarman(0.001, 0.25), fcross; leak=0.999)
        Nrepeat = 10
        Ntimeseries = 2048
        psd_ol_atm, psd_ol_fastncp, psd_ol_slowncp, psd_ol_fastnoise, psd_ol_slownoise = zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints)
        psd_cl_all, psd_cl_atm, psd_cl_fastncp, psd_cl_slowncp, psd_cl_fastnoise, psd_cl_slownoise = zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints)
        etf_atm, etf_fastncp, etf_slowncp, etf_fastnoise, etf_slownoise = zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints)
        @showprogress for _ in 1:Nrepeat
            ts = generate_openloop_timeseries(sim, Ntimeseries)
            slowphase_ts = timedomain_closedloop(sim, (atm=ts.atm, fastncp=ts.fastncp, slowncp=ts.slowncp, fastnoise=ts.fastnoise, slownoise=ts.slownoise/sqrt(sim.R)))
            slowphase_ts_atm = timedomain_closedloop(sim, (atm=ts.atm, fastncp=ts.fastncp*0, slowncp=ts.slowncp*0, fastnoise=ts.fastnoise*0, slownoise=ts.slownoise*0))
            slowphase_ts_fastncp = timedomain_closedloop(sim, (atm=ts.atm*0, fastncp=ts.fastncp, slowncp=ts.slowncp*0, fastnoise=ts.fastnoise*0, slownoise=ts.slownoise*0))
            slowphase_ts_slowncp = timedomain_closedloop(sim, (atm=ts.atm*0, fastncp=ts.fastncp*0, slowncp=ts.slowncp, fastnoise=ts.fastnoise*0, slownoise=ts.slownoise*0))
            slowphase_ts_fastnoise = timedomain_closedloop(sim, (atm=ts.atm*0, fastncp=ts.fastncp*0, slowncp=ts.slowncp*0, fastnoise=ts.fastnoise, slownoise=ts.slownoise*0))
            slowphase_ts_slownoise = timedomain_closedloop(sim, (atm=ts.atm*0, fastncp=ts.fastncp*0, slowncp=ts.slowncp*0, fastnoise=ts.fastnoise*0, slownoise=ts.slownoise/sqrt(sim.R)))

            psd_atm = gen_avg_per_unb(ts.atm, db_Nfreqpoints)[2:k]
            psd_fastncp = gen_avg_per_unb(ts.fastncp, db_Nfreqpoints)[2:k]
            psd_slowncp = gen_avg_per_unb(ts.slowncp, db_Nfreqpoints)[2:k]
            psd_fastnoise = gen_avg_per_unb(ts.fastnoise, db_Nfreqpoints)[2:k]
            psd_slownoise = gen_avg_per_unb(ts.slownoise / sqrt(sim.R), db_Nfreqpoints)[2:k]

            psd_cl_all_i = gen_avg_per_unb(slowphase_ts, db_Nfreqpoints)[2:k]
            psd_cl_atm_i = gen_avg_per_unb(slowphase_ts_atm, db_Nfreqpoints)[2:k]
            psd_cl_fastncp_i = gen_avg_per_unb(slowphase_ts_fastncp, db_Nfreqpoints)[2:k]
            psd_cl_slowncp_i = gen_avg_per_unb(slowphase_ts_slowncp, db_Nfreqpoints)[2:k]
            psd_cl_fastnoise_i = gen_avg_per_unb(slowphase_ts_fastnoise, db_Nfreqpoints)[2:k]
            psd_cl_slownoise_i = gen_avg_per_unb(slowphase_ts_slownoise, db_Nfreqpoints)[2:k]

            psd_ol_atm += psd_atm
            psd_ol_fastncp += psd_fastncp
            psd_ol_slowncp += psd_slowncp
            psd_ol_fastnoise += psd_fastnoise
            psd_ol_slownoise += psd_slownoise

            psd_cl_all += psd_cl_all_i
            psd_cl_atm += psd_cl_atm_i
            psd_cl_fastncp += psd_cl_fastncp_i
            psd_cl_slowncp += psd_cl_slowncp_i
            psd_cl_fastnoise += psd_cl_fastnoise_i
            psd_cl_slownoise += psd_cl_slownoise_i

            etf_atm += psd_cl_atm_i ./ psd_atm
            etf_fastncp += psd_cl_fastncp_i ./ psd_fastncp
            etf_slowncp += psd_cl_slowncp_i ./ psd_slowncp
            etf_fastnoise += psd_cl_fastnoise_i ./ psd_fastnoise
            etf_slownoise += psd_cl_slownoise_i ./ psd_slownoise
        end
        psd_ol_atm /= Nrepeat
        psd_ol_fastncp /= Nrepeat
        psd_ol_slowncp /= Nrepeat
        psd_ol_fastnoise /= Nrepeat
        psd_ol_slownoise /= Nrepeat

        psd_cl_all /= Nrepeat
        psd_cl_atm /= Nrepeat
        psd_cl_fastncp /= Nrepeat
        psd_cl_slowncp /= Nrepeat
        psd_cl_fastnoise /= Nrepeat
        psd_cl_slownoise /= Nrepeat

        etf_atm /= Nrepeat
        etf_fastncp /= Nrepeat
        etf_slowncp /= Nrepeat
        etf_fastnoise /= Nrepeat
        etf_slownoise /= Nrepeat

        c = rand(Int)
        Plots.plot!(sim.fr, error_at_f_X.(sim.fr, Ref(sim)), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="Closed-loop PSD", label="Frequency, fcross=$fcross Hz", legend=:bottomleft, color=c, ls=:dash)
        Plots.plot!(freq, psd_cl_all, label="Time, fcross=$fcross Hz", color=c)
    end
    p
end

begin
    optpars = load("data/all_opt.jld2")
    xerrs = Dict()
    r0_ncps = vcat(0.6:0.2:2.0, 3:6)
    for f_crossover in [50.0, 100.0, 200.0, 500.0]        
        @showprogress for r0_ncp in r0_ncps
            this_xerrs = []
            Nrepeat = 100
            Ntimeseries = 2048
            vk_ncp = VonKarman(0.001, 0.25 * (r0_ncp)^(-5/3))
            sim = simgen_ichpf(optpars["($r0_ncp, $f_crossover)"][3]..., vk_ncp, f_crossover; leak=0.999)
            xerrs["($r0_ncp, $f_crossover) freq"] = notched_error_X(sim)
            for _ in 1:Nrepeat
                ts = generate_openloop_timeseries(sim, Ntimeseries)
                slowphase_ts = timedomain_closedloop(sim, (atm=ts.atm, fastncp=ts.fastncp, slowncp=ts.slowncp, fastnoise=ts.fastnoise, slownoise=ts.slownoise/sqrt(sim.R)))
                psd_cl_all_i = gen_avg_per_unb(slowphase_ts, db_Nfreqpoints)[2:k]
                push!(
                    this_xerrs,
                    sqrt(sum(diff(freq)[1] * psd_cl_all_i[2:end]))
                )
            end
            xerrs["($r0_ncp, $f_crossover)"] = Statistics.quantile(this_xerrs, [0.25, 0.5, 0.75])
        end
    end
    
    hairdryer_td = []
    for f_crossover in [50.0, 100.0, 200.0, 500.0]
        p = Plots.plot(title="fcross = $f_crossover Hz", xlabel="r₀ NCP (m)", ylabel="CL error at X (rad)", xscale=:log10, xticks=([0.6, 1.0, 1.4, 2.0, 3.0, 4.0, 6.0], [0.6, 1.0, 1.4, 2.0, 3.0, 4.0, 6.0]))
        xerrs_pl, errplus, errminus = Float64[], Float64[], Float64[]
        for r0_ncp in r0_ncps
            xerrs_q = xerrs["($r0_ncp, $f_crossover)"]
            push!(xerrs_pl, xerrs_q[2])
            push!(errplus, xerrs_q[3] - xerrs_q[2])
            push!(errminus, xerrs_q[2] - xerrs_q[1])
        end
        Plots.scatter!(r0_ncps, xerrs_pl, yerror=(errminus, errplus), color=1, legend=nothing)
        Plots.plot!(r0_ncps, [xerrs["($r0_ncp, $f_crossover) freq"] for r0_ncp in r0_ncps], color=2, msw=0, legend=nothing)
        push!(hairdryer_td, p)
    end
    Plots.plot(hairdryer_td...)
end

