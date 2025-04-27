using multiwfs
using NPZ
using CairoMakie

sim = simgen_ichpf(0.4, 1.4, 15.0, VonKarman(0.001, 0.25), 50.0; leak=0.999)

db_Nfreqpoints = 1000
Nfreqpoints = db_Nfreqpoints ÷ 2
k = Nfreqpoints + 1
freq = range(0.0, 500.0, k)[2:end]

begin
    Nrepeat = 100
    Ntimeseries = 2048
    psd_ol_atm, psd_ol_fastncp, psd_ol_slowncp, psd_ol_fastnoise, psd_ol_slownoise = zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints)
    psd_cl_atm, psd_cl_fastncp, psd_cl_slowncp, psd_cl_fastnoise, psd_cl_slownoise = zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints)
    etf_atm, etf_fastncp, etf_slowncp, etf_fastnoise, etf_slownoise = zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints)
    for _ in 1:Nrepeat
        ts = generate_openloop_timeseries(sim, Ntimeseries)
        slowphase_ts_atm = timedomain_closedloop(sim, (atm=ts.atm, fastncp=ts.fastncp*0, slowncp=ts.slowncp*0, fastnoise=ts.fastnoise*0, slownoise=ts.slownoise*0))
        slowphase_ts_fastncp = timedomain_closedloop(sim, (atm=ts.atm*0, fastncp=ts.fastncp, slowncp=ts.slowncp*0, fastnoise=ts.fastnoise*0, slownoise=ts.slownoise*0))
        slowphase_ts_slowncp = timedomain_closedloop(sim, (atm=ts.atm*0, fastncp=ts.fastncp*0, slowncp=ts.slowncp, fastnoise=ts.fastnoise*0, slownoise=ts.slownoise*0))
        slowphase_ts_fastnoise = timedomain_closedloop(sim, (atm=ts.atm*0, fastncp=ts.fastncp*0, slowncp=ts.slowncp*0, fastnoise=ts.fastnoise, slownoise=ts.slownoise*0))
        slowphase_ts_slownoise = timedomain_closedloop(sim, (atm=ts.atm*0, fastncp=ts.fastncp*0, slowncp=ts.slowncp*0, fastnoise=ts.fastnoise*0, slownoise=ts.slownoise))

        psd_atm = gen_avg_per_unb(ts.atm, db_Nfreqpoints)[2:k]
        psd_fastncp = gen_avg_per_unb(ts.fastncp, db_Nfreqpoints)[2:k]
        psd_slowncp = gen_avg_per_unb(ts.slowncp, db_Nfreqpoints)[2:k]
        psd_fastnoise = gen_avg_per_unb(ts.fastnoise, db_Nfreqpoints)[2:k]
        psd_slownoise = gen_avg_per_unb(ts.slownoise, db_Nfreqpoints)[2:k]

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
    psd_ol_slownoise /= sim.R # adjust for only taking every 10th sample

    psd_cl_atm /= Nrepeat
    psd_cl_fastncp /= Nrepeat
    psd_cl_slowncp /= Nrepeat
    psd_cl_fastnoise /= Nrepeat
    psd_cl_slownoise /= Nrepeat
    psd_cl_slownoise /= sim.R # adjust for only taking every 10th sample

    etf_atm /= Nrepeat
    etf_fastncp /= Nrepeat
    etf_slowncp /= Nrepeat
    etf_fastnoise /= Nrepeat
    etf_slownoise /= Nrepeat

    psd_ol_atm_theory = psd_von_karman.(sim.fr, Ref(sim.vk_atm))
    psd_ol_ncp_theory = psd_von_karman.(sim.fr, Ref(sim.vk_ncp))
    psd_ol_noise_theory = repeat([psd_von_karman(sim.f_noise_crossover, sim.vk_atm)], length(sim.fr))

    fig_ol = Figure()
    axes = [Axis(fig_ol[r,c], xscale=log10, yscale=log10, xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)") for (r, c) in [(1,1), (1,2), (1,3), (2,2), (2,3)]]
    for (ax, f, psdol_measured, psdol_theory) in zip(
        axes,
        [phi_to_X, Lfast_to_X, Lslow_to_X, Nfast_to_X, Nslow_to_X], 
        [psd_ol_atm, psd_ol_fastncp, psd_ol_slowncp, psd_ol_fastnoise, psd_ol_slownoise],
        [psd_ol_atm_theory, psd_ol_ncp_theory, psd_ol_ncp_theory, psd_ol_noise_theory, psd_ol_noise_theory]
    )
        lines!(ax, freq, psdol_measured)
        lines!(ax, sim.fr, psdol_theory, linestyle=:dash, color=:black, label="$(string(f))")
        #CairoMakie.ylims!(ax, 1e-5, 2e3)
    end
    fig_ol
end

begin
    fig_cl = Figure()
    axes = [Axis(fig_cl[r,c], xscale=log10, yscale=log10, xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)") for (r, c) in [(1,1), (1,2), (1,3), (2,2), (2,3)]]
    for (ax, f, psdcl_measured, psdol_theory) in zip(
        axes, 
        [phi_to_X, Lfast_to_X, Lslow_to_X, Nfast_to_X, Nslow_to_X], 
        [psd_cl_atm, psd_cl_fastncp, psd_cl_slowncp, psd_cl_fastnoise, psd_cl_slownoise],
        [psd_ol_atm_theory, psd_ol_ncp_theory, psd_ol_ncp_theory, psd_ol_noise_theory, psd_ol_noise_theory]
    )
        lines!(ax, freq, psdcl_measured)
        lines!(ax, sim.fr, psdol_theory .* abs2.(f.(sim.sT, Ref(sim))), linestyle=:dash, color=:black, label="$(string(f))")
        #CairoMakie.ylims!(ax, 1e-5, 2e3)
    end
    fig_cl
end

begin
    fig = Figure()
    axes = [Axis(fig[r,c], xscale=log10, yscale=log10, xlabel="Frequency (Hz)", ylabel="ETF") for (r, c) in [(1,1), (1,2), (1,3), (2,2), (2,3)]]
    for (ax, f, etf) in zip(axes, [phi_to_X, Lfast_to_X, Lslow_to_X, Nfast_to_X, Nslow_to_X], [etf_atm, etf_fastncp, etf_slowncp, etf_fastnoise, etf_slownoise])
        lines!(ax, freq, etf)
        lines!(ax, sim.fr, abs2.(f.(sim.sT, Ref(sim))), linestyle=:dash, color=:black, label="$(string(f))")
        CairoMakie.ylims!(ax, 1e-5, 2e3)
        axislegend(ax)
    end
    fig
end