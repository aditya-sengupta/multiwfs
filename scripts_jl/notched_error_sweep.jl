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
sr = 2π .* im .* fr / f_loop
zr = exp.(sr)
no_filter = ZPKFilter(0, 0, 1)
vk_atm = VonKarman(v=10)
vk_ncp = VonKarman(v=0.1, rms_target=0.8)

function notched_errors(gain_slow, gain_fast, f_cutoff)
    ar1_high = ar1_filter(f_cutoff, f_loop, "high")
    sys_high = AOSystem(f_loop, 1.0, gain_fast, 0.999, 10, ar1_high)
    sys_low = AOSystem(f_loop / 10, 0.0, gain_slow, 0.999, 1, no_filter)
    if !is_stable(sys_high) || !is_stable(sys_low)
        return Inf, Inf
    end
    Cfast = z -> Hcont(sys_high, log(z)) * Hfilter(sys_high, log(z))
    Cslow = z -> Hcont(sys_low, log(z))
    notched_error_X(Cfast, Cslow, 10, vk_atm, vk_ncp), notched_error_Y(Cfast, Cslow, 10, vk_atm, vk_ncp)
end


notched_errors(1.4, 0.4, 15)

gains_slow = 0.0:0.1:2.0
gains_fast = 0.0:0.1:1.0
cutoff_freqs = 0.0:2.0:100.0

X_errors = zeros(length(gains_slow), length(gains_fast), length(cutoff_freqs))
Y_errors = zeros(length(gains_slow), length(gains_fast), length(cutoff_freqs))

@showprogress for i in 1:length(gains_slow)
    for j in 1:length(gains_fast)
        @threads for k in 1:length(cutoff_freqs)
            X_errors[i,j,k], Y_errors[i,j,k] = notched_errors(gains_slow[i], gains_fast[j], cutoff_freqs[k])
        end
    end
end

minimum(X_errors)
idx_xmin = argmin(X_errors)
gains_slow[idx_xmin[1]], gains_fast[idx_xmin[2]], cutoff_freqs[idx_xmin[3]]
Y_errors[argmin(X_errors)]

minimum(Y_errors)
idx_ymin = argmin(Y_errors)
gains_slow[idx_ymin[1]], gains_fast[idx_ymin[2]], cutoff_freqs[idx_ymin[3]]
X_errors[argmin(Y_errors)]