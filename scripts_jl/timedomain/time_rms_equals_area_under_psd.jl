using multiwfs
using NPZ
using CairoMakie

rms

atm = npzread("data/olt/olt_atm.npy")
atm = rand(100_000)
freqs, psd = genpsd(atm, 1000.0);
sqrt(sum(psd)), std(atm)

begin
    fig = Figure()
    ax = Axis(fig[1,1], xscale=log10, yscale=log10, xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)")
    lines!(ax, freqs, psd)
    fig
end