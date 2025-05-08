using multiwfs
using CairoMakie

vk_atm = VonKarman(1.0, 0.25 * (0.1031)^(-5/3))
vk_ncp = VonKarman(0.001, 0.25)

freq = exp10.(0.0:0.001:log10(500.0))

begin
    fig = Figure()
    ax = Axis(fig[1,1], xscale=log10, yscale=log10)
    lines!(ax, freq, psd_von_karman.(freq, Ref(vk_atm)))
    lines!(ax, freq, psd_von_karman.(freq, Ref(vk_ncp)))
    fig
end

sum(psd_von_karman.(freq, Ref(vk_atm))) / sum(psd_von_karman.(freq, Ref(vk_ncp)))