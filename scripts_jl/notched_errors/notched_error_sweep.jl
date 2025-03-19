using multiwfs
using multiwfs: Hcont, Hfilter, is_stable
using Plots
using Base: product
using Base.Meta: parse
using Base.Threads: @threads
using ProgressMeter: @showprogress
using Optim
using ForwardDiff

# min |X / phi|^2 |PSD atm| + (|X / Lfast|^2 + |X / Lslow|^2) |PSD NCP| 
# subject to |Y / phi|^2 |PSD atm| + (|Y / Lfast|^2 + |Y / Lslow|^2) |PSD NCP| < threshold

f_loop = 1000.0
fr = 10 .^ (-2:0.01:log10(f_loop/2))
sr = 2π .* im .* fr
sT = sr / f_loop
no_filter = ZPKFilter(0, 0, 1)
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * (0.1031)^(-5/3))
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25)
f_noise_crossover = 200.0
noise_normalization = psd_von_karman(200.0, vk_atm)

function systems(gain_slow, gain_fast, f_cutoff)
    ar1_high = ar1_filter(f_cutoff, f_loop, "high")
    # f_loop, frame_delay, gain, leak, fpf
    sys_fast = AOSystem(f_loop, 1.0, gain_fast, 0.999, 1, ar1_high)
    sys_slow = AOSystem(f_loop / 10, 0.1, gain_slow, 0.999, 1, no_filter)
    return sys_fast, sys_slow
end

function notched_errors(gain_slow, gain_fast, f_cutoff)
    λ = 1e3
    penalty = 0
    if gain_slow < 0 || gain_fast < 0 || f_cutoff < 0 || gain_slow > 2.0
        penalty = 1
    end
    sys_fast, sys_slow = systems(gain_slow, gain_fast, f_cutoff)
    Cfast = sT -> Hcont(sys_fast, sT * f_loop) * Hfilter(sys_fast, sT * f_loop)
    Cslow = sT -> Hcont(sys_slow, sT * f_loop)
    if !is_stable(Vector{AOSystem}([sys_fast, sys_slow]))
       penalty = 1
    end
    err_X, err_Y = notched_error_X(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover, f_min=1.0), notched_error_Y(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover, f_min=1.0)
    return err_X + λ * penalty, err_Y + λ * penalty
end


notched_errors(1.4, 0.4, 15.0)

notched_errors(0.0, 0.42, 1.4)
notched_errors(2.0, 0.53, 34.295)
notched_errors(0.1, 0.3, 80.0)

margins(Vector{AOSystem}([systems(3.014, 0.531, 38.246)...]))

gains_slow = 0.0:0.01:2.0
gains_fast = 0.0:0.01:1.0
cutoff_freqs = 0.0:1.0:20.0

"""X_errors = zeros(length(gains_slow), length(gains_fast), length(cutoff_freqs))
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
minimum(X_errors)"""

systems_lisa = systems(1.4, 0.4, 15.0)
systems_opt = systems(2.0, 0.45, 20.0)

res = optimize(x -> notched_errors(x...)[1], [1.4, 0.4, 15.0], NelderMead())
res_params = copy(res.minimizer)
systems_opt = systems(res_params...)

nyquist_plot(Vector{AOSystem}([systems_opt...]))

notched_errors(1.4, 0.4, 15.0)
notched_errors(res_params...)

sTs = 2π .* im * fr / f_loop
begin
    p = plot(xscale=:log10, yscale=:log10, legend=:bottomright, ylims=(1e-7, 1e1), xlabel="Frequency (Hz)", ylabel="Power")
    hline!([1], color=:black, label=nothing)
    for (systems, lstyle) in zip([systems_lisa, systems_opt],[:solid, :dash])
        Cfast = sT -> Hcont(systems[1], sT * f_loop) * Hfilter(systems[1], sT * f_loop)
        Cslow = sT -> Hcont(systems[2], sT * f_loop)
        plot!(fr, abs2.(Lfast_to_X.(sT, Cfast, Cslow, 10)), ls=lstyle, label="xtf_ncp0")
        plot!(fr, abs2.(Lslow_to_X.(sT, Cfast, Cslow, 10)), ls=lstyle, label="xtf_ncp1")
    end
    p
end