using multiwfs
using Plots

atm_aosim = VonKarman(0.3 * 10.0 / 3.0, 0.25 * (0.1031)^(-5/3))
ncp_aosim = VonKarman(0.3 * 0.01 / 3.0, 0.25)

fr = exp10.(-2:0.01:log10(500))

begin
    plot(fr, psd_von_karman.(fr, Ref(atm_aosim)), xscale=:log10, yscale=:log10, yticks=[1e1, 1e-1, 1e-3, 1e-5, 1e-7], ylims=(1e-7, 1e2), label="atm", xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)")
    plot!(fr, psd_von_karman.(fr, Ref(ncp_aosim)), xscale=:log10, yscale=:log10, label="ncp")
    hline!([psd_von_karman(200, atm_aosim)], label="noise")
end