using multiwfs
using multiwfs: Hlint, Hcont, Hfilter
using SciPy
using Plots

# HPF critical frequency of 15 Hz; Chebyshev cutoff factor of 1/20, alpha = 0.91
# Slow sensor gain = 1.4 (runs 10x slower with the same computational delay), fast sensor gain = 0.4

f_loop = 1000.0
f_cutoff = 15.0
ch_rp = -20*log10(0.91)
ch_omegan = [f_cutoff / (f_loop / 2)] # critical frequencies, nyquist = 1
lisa_hpf = ZPKFilter(signal.cheby1(2, ch_rp, ch_omegan, "highpass", output="zpk")...)
no_filter = ZPKFilter(0, 0, 1)

sys_high = AOSystem(f_loop, 1, 0.4, 1.0, 1, lisa_hpf)
sys_low = AOSystem(f_loop, 1, 1.4, 1.0, 10, no_filter)

fr = exp10.(-2:0.01:log10(f_loop/2))
sr = 2π .* im .* fr ./ f_loop
zr = exp.(sr)

plot(fr, abs2.(Hrej.(Ref(sys_high), fr)), xscale=:log10, yscale=:log10)

# moving this over from the LQG-first design notebook
begin
    ar1_low = ar1_filter(f_cutoff, f_loop / 10, "low")
    sys_low = AOSystem(f_loop, 1.0, 0.1, 0.9999999, 10, ar1_low)
    Cslow = z -> (real(log(z) / (2π * im)) < 1/20) ? (Hfilter(sys_low, log(z)) * Hcont(sys_low, log(z))) : 1.0
    plot(fr, abs2.(Lslow_to_Y.(zr, Cfast, Cslow, 10)), xscale=:log10, yscale=:log10)
end
