using multiwfs
using SciPy
using Plots

f_loop = 200.0
f_cutoff = 3.0
ar1_high = ar1_filter(f_cutoff, f_loop, "high")
ch_rp = -20 * log10(0.95)
ch_omegan = [f_cutoff / (f_loop / 2)] # critical frequencies, nyquist = 1
@assert 0 < ch_omegan[1] < 1
cheb2_high = ZPKFilter(signal.cheby1(2, ch_rp, ch_omegan, "highpass", output="zpk")...)
cheb4_high = ZPKFilter(signal.cheby1(4, ch_rp, ch_omegan, "highpass", output="zpk")...)

begin
    fr = (6e-4:1e-4:0.5) .* f_loop
    sr = f2s.(fr ./ f_loop)
    plot(fr, abs2.(transfer_function.(Ref(ar1_high), sr)), xscale=:log10, yscale=:log10, ylim=(1e-5, 1.1), xlabel="Frequency (Hz)", ylabel="Power (normalized)", label="AR(1) HPF", xticks=[1e-1, 1e0, 1e1, 1e2])
    plot!(fr, abs2.(transfer_function.(Ref(cheb2_high), sr)), xscale=:log10, yscale=:log10, label="Chebyshev order 2 HPF")
    plot!(fr, abs2.(transfer_function.(Ref(cheb4_high), sr)), xscale=:log10, yscale=:log10, label="Chebyshev order 4 HPF")
end