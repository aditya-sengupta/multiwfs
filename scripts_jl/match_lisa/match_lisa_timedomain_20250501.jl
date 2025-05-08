using multiwfs
using NPZ
using CairoMakie

sim_lisa = simgen_ichpf(0.4, 1.4, 15.0, VonKarman(0.001, 0.25 * (0.6)^(-5/3)), 500.0; leak=0.999)

db_Nfreqpoints = 1000
Nfreqpoints = db_Nfreqpoints ÷ 2
k = Nfreqpoints + 1
freq = range(0.0, 500.0, k)[2:end]

Ntimeseries = 20480
ts = generate_openloop_timeseries(sim_lisa, Ntimeseries, [1, 0, 0, 0, 0])
slowphase_ts = timedomain_closedloop(sim_lisa, ts)

begin
    fig = Figure()
    ax = Axis(fig[1,1], xscale=log10, yscale=log10, xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)")
    lines!(ax, freq, gen_avg_per_unb(ts.atm, 500) .* abs2.(phi_to_X.(2π * im * freq / sim_lisa.f_loop, Ref(sim_lisa)))
    , label="Time-domain OL PSD times analytic ETF")
    lines!(ax, freq, gen_avg_per_unb(slowphase_ts, 500), label="Time-domain CL PSD")
    axislegend(ax)
    fig
end

begin
    fig = Figure()
    ax = Axis(fig[1,1], xscale=log10, yscale=log10, xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)")
    lines!(ax, freq, abs2.(phi_to_X.(2π * im * freq / sim_lisa.f_loop, Ref(sim_lisa)))
    , label="Analytic ETF")
    lines!(ax, freq, gen_avg_per_unb(slowphase_ts, 500) ./ gen_avg_per_unb(ts.atm, 500), label="Time-domain ETF")
    axislegend(ax, position=:rb)
    fig
end