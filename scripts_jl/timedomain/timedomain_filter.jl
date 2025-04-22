using multiwfs
using CairoMakie
using StatsBase: mean

filt = ar1_filter(15.0, 1000.0, "high")
N = 100_000
x = cumsum(2 * rand(N) .- 1)
y = zeros(N)
reset!(filt)
for i in 1:N
    y[i] = output!(filt, x[i])[1]
end
freqs, psd_unfiltered = genpsd(x, 1000.0)
freqs, psd_filtered = genpsd(y, 1000.0)

begin
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Frequency (Hz)", ylabel="PSD", xscale=log10, yscale=log10, title="Frequency response of filter")
    lines!(ax, freqs, psd_unfiltered, psd_filtered, label="Unfiltered")
    lines!(ax, freqs, psd_filtered, psd_filtered, label="Filtered")
    axislegend(ax)
    fig
end

begin
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Frequency (Hz)", ylabel="PSD", xscale=log10, yscale=log10, title="Frequency response of filter")
    lines!(ax, freqs, psd_filtered ./ psd_unfiltered, label="Time-domain ETF")
    lines!(ax, freqs, abs2.(transfer_function.(Ref(filt), 2π * im * freqs / 1000.0)), label="Analytic ETF")
    axislegend(ax)
    fig
end