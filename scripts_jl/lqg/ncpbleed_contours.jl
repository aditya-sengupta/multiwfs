using multiwfs
using Plots

f_loop = 1000.0
fr = exp10.(-2:0.01:log10(f_loop/2))
sr = 2π * im * fr ./ f_loop
zr = exp.(sr)

Cabs = exp10.(-4:0.01:0)
Cphi = 0:(π/200):2π

f_cutoff = 0.05
ar1_low = ar1_filter(f_cutoff, f_loop / 10, "low")
sys_low_regular = AOSystem(f_loop, 1.0, 0.1, 1 - 1e-7, 10, ar1_low)
search_gain!(sys_low_regular)
sys_low_rejectlslow = AOSystem(f_loop, 1.0, 0.001, 1 - 1e-2, 10, ar1_low)
search_gain!(sys_low_rejectlslow)
sys_low_rejectlslow.gain /= 5
Cslow_regular = z -> (real(log(z) / (2π * im)) < 1/20) ? (Hfilter(sys_low_regular, log(z)) * Hcont(sys_low_regular, log(z))) : 1.0

Cslow_rejectlslow = z -> (real(log(z) / (2π * im)) < 1/20) ? (Hfilter(sys_low_rejectlslow, log(z)) * Hcont(sys_low_rejectlslow, log(z))) : 1.0

begin
    Lslow_to_Y_abs = abs.(hcat([[Lslow_to_Y(zr[1], z -> Cabs_v * exp(1im * Cphi_v), Cslow_regular, 10) for Cabs_v in Cabs] for Cphi_v in Cphi]...))
    heatmap(Cabs, Cphi, Lslow_to_Y_abs, xscale=:log10, title="Min slow NCP -> fast WFS rejection = $(round(minimum(Lslow_to_Y_abs), digits=5))", xlabel="|Cfast|", ylabel="arg Cfast (rad)")
end

begin
    Lslow_to_Y_abs = abs.(hcat([[Lslow_to_Y(zr[1], z -> Cabs_v * exp(1im * Cphi_v), Cslow_rejectlslow, 10) for Cabs_v in Cabs] for Cphi_v in Cphi]...))
    heatmap(Cabs, Cphi, Lslow_to_Y_abs, xscale=:log10, title="Min slow NCP -> fast WFS rejection = $(round(minimum(Lslow_to_Y_abs), digits=5))", xlabel="|Cfast|", ylabel="arg Cfast (rad)")
end


begin
    Cfast_const = 1e-1
    plot(fr, abs2.(Lslow_to_X.(zr, z -> Cfast_const, Cslow, 10)), xscale=:log10, xlabel="Frequency (Hz)", ylabel="|NCP-bleed TFs|²", yscale=:log10, label="Lslow -> X, don't reject Lslow -> Y", legend=:left, ylim=(1e-8, 10), yticks=[1e-8, 1e-6, 1e-4, 1e-2, 1e0], color=1)
    plot!(fr, abs2.(Lslow_to_Y.(zr, z -> Cfast_const, Cslow, 10)), xscale=:log10, xlabel="Frequency (Hz)", ylabel="|NCP-bleed TFs|²", yscale=:log10, label="Lslow -> Y, don't reject Lslow -> Y", color=1, ls=:dash)
    plot!(fr, abs2.(Lslow_to_X.(zr, z -> Cfast_const, Cslow_rejectlslow, 10)), xscale=:log10, label="Lslow -> X, do reject Lslow -> Y", color=2)
    plot!(fr, abs2.(Lslow_to_Y.(zr, z -> Cfast_const, Cslow_rejectlslow, 10)), xscale=:log10, label="Lslow -> Y, do reject Lslow -> Y", color=2, ls=:dash)
end
