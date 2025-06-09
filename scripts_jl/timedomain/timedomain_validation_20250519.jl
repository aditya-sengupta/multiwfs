using multiwfs
using NPZ
using CairoMakie

# if we do moving average instead of sample and hold
# does the mismatch go away?

begin
    sim = simgen_ichpf(0.4, 1.4, 15.0, VonKarman(0.001, 0.25), 50.0; leak=0.999)
    db_Nfreqpoints = 1000
    Nfreqpoints = db_Nfreqpoints ÷ 2
    k = Nfreqpoints + 1
    freq = range(0.0, 500.0, k)[2:end]
    slow_nyquist = sim.f_loop / (2 * sim.R)
    slownoise_modifier = freq .<= slow_nyquist
    activations = [1, 1, 1, 1, 1]
    # activations = [0, 0, 0, 0, 1]
    Ntimeseries = 204800
    psd_ol_atm = zeros(Nfreqpoints)
    psd_cl_atm = zeros(Nfreqpoints)
    etf_atm = zeros(Nfreqpoints)
    ts = generate_openloop_timeseries(sim, Ntimeseries)
    slowphase_ts = timedomain_closedloop(sim, (atm=ts.atm*activations[1], fastncp=ts.fastncp*activations[2], slowncp=ts.slowncp*activations[3], fastnoise=ts.fastnoise*activations[4], slownoise=ts.slownoise*activations[5]/sqrt(sim.R)))

    # time-domain curves for OL and CL PSDs
    psd_ol_atm = gen_avg_per_unb(ts.atm, db_Nfreqpoints)[2:k]
    psd_ol_fastncp = gen_avg_per_unb(ts.fastncp, db_Nfreqpoints)[2:k]
    psd_ol_slowncp = gen_avg_per_unb(ts.slowncp, db_Nfreqpoints)[2:k]
    psd_ol_fastnoise = gen_avg_per_unb(ts.fastnoise, db_Nfreqpoints)[2:k]
    psd_ol_slownoise = gen_avg_per_unb(ts.slownoise / sqrt(sim.R), db_Nfreqpoints)[2:k]
    psd_ol_measured = psd_ol_atm * activations[1] + psd_ol_fastncp * activations[2] + psd_ol_slowncp * + activations[3] + psd_ol_fastnoise * activations[4] + psd_ol_slownoise * activations[5] .* slownoise_modifier
    psd_cl_measured = gen_avg_per_unb(slowphase_ts, db_Nfreqpoints)[2:k]

    # ideal curves for OL and CL PSDs
    sT = 2π * im * freq / sim.f_loop
    psd_ol_atm_theory = psd_von_karman.(sim.fr, Ref(sim.vk_atm))
    psd_ol_ncp_theory = psd_von_karman.(sim.fr, Ref(sim.vk_ncp))
    psd_ol_noise_theory = repeat([psd_von_karman(sim.f_noise_crossover, sim.vk_atm)], length(sim.fr))
    psd_ol_theory = psd_ol_atm * activations[1] + psd_ol_fastncp * activations[2] + psd_ol_slowncp * activations[3] + psd_ol_fastnoise * activations[4] + psd_ol_slownoise * activations[5]
    psd_cl_theory = psd_ol_atm .* abs2.(phi_to_X.(sT, Ref(sim))) * activations[1] + psd_ol_fastncp .* abs2.(Lfast_to_X.(sT, Ref(sim))) * activations[2] + psd_ol_slowncp .* abs2.(Lslow_to_X.(sT, Ref(sim))) * activations[3] + psd_ol_fastnoise .* abs2.(Nfast_to_X.(sT, Ref(sim))) * activations[4] + psd_ol_slownoise .* abs2.(Nslow_to_X.(sT, Ref(sim))) / sim.R * activations[5]

    fig_olcl = Figure()
    ax = Axis(fig_olcl[1,1], xscale=log10, yscale=log10, xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)", title="(atm, fast NCP, slow NCP, fast noise, slow noise) = $activations")
    lines!(ax, freq, psd_ol_measured, color=:teal, label="Open loop")
    lines!(ax, freq, psd_cl_measured, color=:orange, label="Closed loop")
    lines!(ax, freq, psd_cl_theory, linestyle=:dash, color=:orange)
    vlines!(ax, [50.0], color=:black, linestyle=:dash)
    axislegend(ax, position=:lb)
    CairoMakie.ylims!(ax, 1e-6, 1e2)
    fig_olcl
end

sqrt(sum(diff(freq)[1] * psd_cl_theory[2:end]))
sqrt(sum(diff(freq)[1] * psd_cl_measured[2:end]))

# PSD of output time-series in CL = each TF * PSD of each input