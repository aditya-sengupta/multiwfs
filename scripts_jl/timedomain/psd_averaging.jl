using multiwfs
using CairoMakie

vk_ncp = VonKarman(0.001, 0.25)
freqs = range(0.0, 500.0, 1025)[2:end]
begin
    periodograms = []
    for _ in 1:100
        ts = generate_time_series(f -> psd_von_karman(f, vk_ncp), 500.0, 2048)
        push!(periodograms, gen_avg_per_unb(ts, 2048)[2:1025])
    end
    Navg = Observable(100)
    periodograms_scaled = periodograms
    avg_periodogram = @lift(sum(periodograms_scaled[1:$Navg]) / $Navg)
    fig = Figure()
    ax = Axis(fig[1,1], xscale=log10, yscale=log10, xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)", title=@lift("Averaging over $($Navg) frames"))
    lines!(ax, freqs, psd_von_karman.(freqs, Ref(vk_ncp)), color=:black, linestyle=:dash)
    lines!(ax, freqs, avg_periodogram)
    fig
end

record(fig, "averaging_periodogram.mp4",vcat(1:5, 10:5:100); framerate=10) do n
    Navg[] = n
end