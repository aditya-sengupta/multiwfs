using multiwfs
using Plots
using Distributions
using SciPy
using multiwfs: ar1_high_filter, ar1_low_filter

# assess Ben's 1 rad rms criterion for AR(1), Cheb2, Cheb4

begin
    f_cutoff = 3.0
    f_loop = 200.0
    ar1_high = ar1_high_filter(f_cutoff, f_loop)
    # parameters are:
    # f_loop, frame_delay, gain, leak, fpf
    sys_high = AOSystem(f_loop, 1.0, 0.15, 0.999, 10, ar1_high[1], ar1_high[2])
    sys_test = AOSystem(f_loop, 0.0, 0.3, 0.999, 10, [1], [1])
    sys_slow = AOSystem(f_loop / 10, 0.0, 1.0, 0.999, 1, [1], [1])
end;


ch_order = 1
 # order of the filter
cheby_cutoff_fac = f_cutoff / (sys_high.f_loop / 2)
ch_rp = -20 * log10(0.95)
@assert 0 < cheby_cutoff_fac < 1
ch_omegan = [cheby_cutoff_fac] # critical frequencies, nyquist = 1
hp_ch_numer, hp_ch_denom = signal.cheby1(ch_order, ch_rp, ch_omegan, "highpass", output="ba")
w, h = signal.freqs(hp_ch_numer, hp_ch_denom)
plot(w, 20 .* log10.(abs.(h)), xscale=:log10, xlabel="Frequency (Hz)", xticks=[1e-2, 1e-1, 1e0, 1e1, 1e2], ylabel="20 log10(power) (dB)")
# order 2 should behave similarly to the AR
# combined TF should have a saddle in both 2 and 4
# order 4 should have no improvement
# in the zpk form, each individual filter term should remember its last state
# in the ba form, you need a length n history

# let's make some open-loop turbulence
begin
    N = 50000
    f_loop = 200
    open_loop_t = zeros(N)
    for i in 2:N
        open_loop_t[i] = 0.995 * open_loop_t[i-1] + rand(Normal(0, 0.001))
    end
    const open_loop = open_loop_t
end

ol_psd_p = psd(open_loop, f_loop)
f, ol_psd = freq(ol_psd_p)[2:end], power(ol_psd_p)[2:end]
etf_regular = power(psd(integrator_control(open_loop, 0.3, 0.999, 1, delay_frames=0)))[2:end] ./ ol_psd
etf_slow = power(psd(integrator_control(open_loop, 1.0, 0.999, 10, delay_frames=1)))[2:end] ./ ol_psd
etf_filt = power(psd(integrator_control(open_loop, 1.0, 0.999, 10, hpf_gain=0.15, delay_frames=1)))[2:end] ./ ol_psd

begin
    plot(title="Estimated ETFs - Slow WFS @ 20 Hz", legend=:bottomright, xticks=[1e-1, 1e0, 1e1, 1e2])
    plot_psd_p!(f, etf_regular, label="Regular g = 0.30", color=:green)
    plot_psd_p!(f, etf_slow, label="Slow gₛ = 1.00", color=2)
    plot_psd_p!(f, etf_filt, label="Slow + Fast-HPF gₛ = 1.00, g = 0.15", color=1)
    plot_psd_p!(f, abs2.(1 ./ (1 .+ Hol_unfiltered.(Ref(sys_test), f))), color=:black,  label=nothing)
    plot_psd_p!(f, abs2.(1 ./ (1 .+ Hol_unfiltered.(Ref(sys_slow), f))), color=:black,  label=nothing)
    plot_psd_p!(f, abs2.(1 ./ (1 .+ Hol_unfiltered.(Ref(sys_slow), f) .+ Hol.(Ref(sys_high), f))), color=:black,  label=nothing)
end