using multiwfs
using multiwfs: Hcont, Hfilter, is_stable
using Plots
using Base: product
using Base.Meta: parse
using Base.Threads: @threads
using ProgressMeter: @showprogress

# min |X / phi|^2 |PSD atm| + (|X / Lfast|^2 + |X / Lslow|^2) |PSD NCP| 
# subject to |Y / phi|^2 |PSD atm| + (|Y / Lfast|^2 + |Y / Lslow|^2) |PSD NCP| < threshold

f_loop = 1000.0
fr = 10 .^ (-2:0.01:log10(f_loop/2))
sr = 2π .* im .* fr
sT = sr / f_loop
no_filter = ZPKFilter(0, 0, 1)
vk_atm = VonKarman(v=10)
vk_ncp = VonKarman(v=0.1, rms_target=0.8)
f_noise_crossover = 200.0
noise_normalization = psd_von_karman(200.0, vk_atm)

function notched_errors(gain_slow, gain_fast, f_cutoff)
    ar1_high = ar1_filter(f_cutoff, f_loop, "high")
    # f_loop, frame_delay, gain, leak, fpf
    sys_fast = AOSystem(f_loop, 1.0, gain_fast, 0.999, 1, ar1_high)
    sys_slow = AOSystem(f_loop / 10, 0.1, gain_slow, 0.999, 1, no_filter)
    Cfast = sT -> Hcont(sys_fast, sT * f_loop) * Hfilter(sys_fast, sT * f_loop)
    Cslow = sT -> Hcont(sys_slow, sT * f_loop)
    if !is_stable(sys_slow)
       return Inf, Inf 
    end
    notched_error_X(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover), notched_error_Y(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover)
end

notched_errors(1.4, 0.4, 15)
notched_errors(1.73, 0.75, 2)
notched_errors(2.0, 0.75, 2.0)

gains_slow = 1.6:0.01:2.0
gains_fast = 0.0:0.05:1.0
cutoff_freqs = 0.0:1.0:5.0

X_errors = zeros(length(gains_slow), length(gains_fast), length(cutoff_freqs))
Y_errors = zeros(length(gains_slow), length(gains_fast), length(cutoff_freqs))

@showprogress for i in eachindex(gains_slow)
    for j in eachindex(gains_fast)
        @threads for k in eachindex(cutoff_freqs)
            X_errors[i,j,k], Y_errors[i,j,k] = notched_errors(gains_slow[i], gains_fast[j], cutoff_freqs[k])
        end
    end
end

minimum(X_errors)
idx_xmin = argmin(X_errors)
gains_slow[idx_xmin[1]], gains_fast[idx_xmin[2]], cutoff_freqs[idx_xmin[3]]
Y_errors[argmin(X_errors)]

X_errors[Y_errors .>= 1.0] .= Inf
argmin(X_errors)
minimum(X_errors)