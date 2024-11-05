using FFTW
using Plots

begin
    dt = 0.001
    t = 0.0:dt:100.0 # s
    N = length(t)
    freq = 1 # Hz
    y = sin.(2π*freq*t)

    y_psd = abs2.(fft(y))
    ω₀ = 1 / (N * dt)
    freqs = fftfreq(length(y), 1/dt)

    plot(
        plot(t[1:1000], y[1:1000], xlabel="Time (s)", ylabel="Signal"),
        plot(freqs[2:length(y)÷2], y_psd[2:length(y)÷2], xscale=:log10, yscale=:log10,
        xlabel="Frequency (Hz)", ylabel="Power")
    )
end