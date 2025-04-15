using FFTW
using DSP
using Plots

# generate PSD for AR process with coefficients Φ
function ar_psd(f, Φ::Vector{Float64})
    order = length(Φ)
    return 1 / abs2(1 - sum(Φ .* exp.(-im * 2π * f * (1:order))))
end

"""
Generate a time series with a given power spectral density (PSD) model.

# Arguments
- `psd_model::Function`: A function that takes a frequency and returns the PSD value at that frequency.
- `fmax::Float64`: The maximum frequency of interest.
- `npoints::Int`: The number of points in the discrete PSD.

# Returns
- `time_series::Vector{Float64}`: The generated time series.

copilot-generated off this:
https://dsp.stackexchange.com/questions/76660/generate-a-time-series-from-power-spectral-density 
"""
function generate_time_series(psd_model::Function, fmax::Float64, npoints::Int)
    fs = 2 * fmax
    df = fs / npoints
    freqs = df * (0:(npoints ÷ 2))
    psd_values = psd_model.(freqs)
    amplitudes = sqrt.(npoints * fs * psd_values)
    phases = 2π * rand(length(amplitudes))
    complex_amplitudes = amplitudes .* exp.(-im .* phases)
    full_spectrum = vcat(complex_amplitudes, conj(reverse(complex_amplitudes[2:end-1])))
    time_series_complex = ifft(full_spectrum)
    time_series = real(time_series_complex)
    return time_series
end

vk_atm = VonKarman(10.0, 0.25 * 0.1031^(-5/3))
vk_ncp = VonKarman(0.01, 0.25 * 0.6^(-5/3))

function compare_psds(vk; order=1)
    ts = generate_time_series(f -> psd_von_karman(f, vk), 500.0, Int(1e6))
    fr = freq(psd(ts, 1000.0))
    X = hcat([ts[i:end-(order+1-i)] for i in 1:order]...)
    y = ts[order+1:end]
    Φ = X \ y
    c = rand(Int)
    ar_psd_values = ar_psd.(fr[2:end] ./ 1000.0, Ref(Φ))
    plot!(fr[2:end], ar_psd_values ./ ar_psd(fr[2] / 1000, Φ), label="AR($order) PSD", yscale=:log10, xscale=:log10, color=c)
    plot!(fr[2:end], psd_von_karman.(fr[2:end], Ref(vk)) / psd_von_karman(fr[2], vk), xscale=:log10, yscale=:log10, label="von Karman PSD", legend=:bottomleft, color=c, ls=:dash)
end

begin
    pl = plot(xlabel="Frequency (Hz)", ylabel="PSD (normalized)")
    for vk in [vk_atm, vk_ncp]
        compare_psds(vk)
    end
    pl
end