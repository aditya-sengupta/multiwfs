using multiwfs
using Plots
pgfplotsx()

vk_ncp = VonKarman(0.001, 0.25)
f_crossover = 50.0
sim = simgen_ichpf(0.4, 1.4, 15.0, vk_ncp, f_crossover)

db_Nfreqpoints = 1000
Nfreqpoints = db_Nfreqpoints ÷ 2
k = Nfreqpoints + 1
freq_td = range(0.0, 500.0, k)[2:end]

begin
    Nrepeat = 1
    Ntimeseries = 2048
    psd_ol_atm, psd_ol_fastncp, psd_ol_slowncp, psd_ol_fastnoise, psd_ol_slownoise = zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints), zeros(Nfreqpoints)
    for _ in 1:Nrepeat
        ts = generate_openloop_timeseries(sim, Ntimeseries)

        psd_atm = gen_avg_per_unb(ts.atm, db_Nfreqpoints)[2:k]
        psd_fastncp = gen_avg_per_unb(ts.fastncp, db_Nfreqpoints)[2:k]
        psd_fastnoise = gen_avg_per_unb(ts.fastnoise, db_Nfreqpoints)[2:k]

        psd_ol_atm += psd_atm
        psd_ol_fastncp += psd_fastncp
        psd_ol_fastnoise += psd_fastnoise
    end
    psd_ol_atm /= Nrepeat
    psd_ol_fastncp /= Nrepeat
    psd_ol_fastnoise /= Nrepeat

    p = plot(sim.fr, psd_von_karman.(sim.fr, Ref(sim.vk_atm)), xscale=:log10, yscale=:log10, xticks=exp10.(-3:2), xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)", label="Atmosphere", ls=:dash, color=1, legend=:topright)
    plot!(freq_td, psd_ol_atm, alpha=0.7, color=1, label="Atmosphere time-domain")
    plot!(sim.fr, psd_von_karman.(sim.fr, Ref(sim.vk_ncp)), label="NCP", ls=:dash, color=2)
    plot!(freq_td, psd_ol_fastncp, alpha=0.7, color=2, label="NCP time-domain")
    plot!(sim.fr, psd_von_karman.(repeat([f_crossover], length(sim.fr)), Ref(sim.vk_atm)), label="Noise", ls=:dash, color=3)
    plot!(freq_td, psd_ol_fastnoise, alpha=0.7, color=3, label="Noise time-domain")

    Plots.savefig(p, "paper/figures/openloop_psds.pdf")
end