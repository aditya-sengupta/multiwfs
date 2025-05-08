using multiwfs
using NPZ
using CairoMakie

sim = simgen_ichpf(0.4, 1.4, 15.0, VonKarman(0.001, 0.25), 100.0; leak=0.999)

db_Nfreqpoints = 100
Nfreqpoints = db_Nfreqpoints ÷ 2
k = Nfreqpoints + 1
freq = range(0.0, 500.0, k)[2:end]

begin
    activations = [1, 0, 0, 0, 0]
    Nrepeat = 1
    Ntimeseries = 2048
    psd_ol = zeros(Nfreqpoints)
    psd_cl = zeros(Nfreqpoints)
    etf = zeros(Nfreqpoints)
    for _ in 1:Nrepeat
        ts = generate_openloop_timeseries(sim, Ntimeseries)
        ts = (atm=ts.atm*activations[1], fastncp=ts.fastncp*activations[2], slowncp=ts.slowncp*activations[3], fastnoise=ts.fastnoise*activations[4], slownoise=ts.slownoise*activations[5])
        slowphase_ts = timedomain_closedloop(sim, ts)

        psd_ol_i = gen_avg_per_unb(ts.atm, db_Nfreqpoints)[2:k] 
        + gen_avg_per_unb(ts.fastncp, db_Nfreqpoints)[2:k] 
        + gen_avg_per_unb(ts.slowncp, db_Nfreqpoints)[2:k] 
        + gen_avg_per_unb(ts.fastnoise, db_Nfreqpoints)[2:k] 
        + gen_avg_per_unb(ts.slownoise, db_Nfreqpoints)[2:k]

        psd_cl_i = gen_avg_per_unb(slowphase_ts, db_Nfreqpoints)[2:k]

        psd_ol += psd_ol_i
        psd_cl += psd_cl_i

        etf += psd_cl_i ./ psd_ol_i
    end
    psd_ol /= Nrepeat
    psd_cl /= Nrepeat
    etf /= Nrepeat

    psd_ol_atm_theory = psd_von_karman.(sim.fr, Ref(sim.vk_atm))
    psd_ol_ncp_theory = psd_von_karman.(sim.fr, Ref(sim.vk_ncp))
    psd_ol_noise_theory = repeat([psd_von_karman(sim.f_noise_crossover, sim.vk_atm)], length(sim.fr))

    ol_theory =
    psd_ol_atm_theory * activations[1]
    + psd_ol_ncp_theory * activations[2]
    + psd_ol_ncp_theory * activations[3]
    + psd_ol_noise_theory * activations[4]
    + psd_ol_noise_theory * activations[5]

    cl_theory =
        psd_ol_atm_theory .* abs2.(phi_to_X.(sim.sT, Ref(sim))) * activations[1]
        + psd_ol_ncp_theory .* abs2.(Lfast_to_X.(sim.sT, Ref(sim))) * activations[2]
        + psd_ol_ncp_theory .* abs2.(Lslow_to_X.(sim.sT, Ref(sim))) * activations[3]
        + psd_ol_noise_theory .* abs2.(Nfast_to_X.(sim.sT, Ref(sim))) * activations[4]
        + psd_ol_noise_theory .* abs2.(Nslow_to_X.(sim.sT, Ref(sim))) * activations[5]

    fig_ol = Figure()
    ax = Axis(fig_ol[1,1], xscale=log10, yscale=log10, xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)")
    lines!(ax, freq, psd_ol)
    lines!(ax, sim.fr, ol_theory, linestyle=:dash, color=:black)
    fig_ol
end

begin
    fig_cl = Figure()
    ax = Axis(fig_cl[1,1], xscale=log10, yscale=log10, xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)")
    lines!(ax, freq, psd_cl)
    lines!(ax, sim.fr, cl_theory, linestyle=:dash, color=:black)
    fig_cl
end

begin
    fig_ol = Figure()
    ax = Axis(fig_ol[1,1], xscale=log10, yscale=log10, xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)")
    lines!(ax, freq, psd_ol)
    lines!(ax, sim.fr, ol_theory, linestyle=:dash, color=:black)
    fig_ol
end