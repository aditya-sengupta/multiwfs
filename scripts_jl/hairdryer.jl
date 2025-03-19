using multiwfs
using ProgressMeter
using Base.Threads: @threads
using Optim
using Plots

r0 = 0.1031
r0_ncp = 1.0
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))

# fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover
f_loop = 1000.0
R = 10

gains_slow = 1.5:0.05:2.0
gains_fast = 0.0:0.05:1.0
cutoff_freqs = 0.0:2.0:80.0

r0_ncp_vals = [0.6, 0.8, 1.0, 1.4, 1.6, 2.0, 4.0, 6.0]

no_filter = ZPKFilter(0, 0, 1)

begin
    opt_xerrs, opt_xerrs_nofilter, opt_xparams = [], [], []
    for r0_ncp in r0_ncp_vals
        vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
        X_errors_nohpf = zeros(length(gains_slow), length(gains_fast))
        X_errors = zeros(length(gains_slow), length(gains_fast), length(cutoff_freqs))
        @showprogress for i in eachindex(gains_slow)
            slow_controller = FilteredIntegrator(gains_slow[i], 0.999, no_filter, R/f_loop)
            for j in eachindex(gains_fast)
                fast_unfiltered = FilteredIntegrator(gains_fast[j], 0.999, no_filter, 1/f_loop)
                sim_unfiltered = Simulation(f_loop, fast_unfiltered, slow_controller, R, vk_atm, vk_ncp, 500.0)
                X_errors_nohpf[i,j] = notched_error_X(sim_unfiltered)
                @threads for k in eachindex(cutoff_freqs)
                    fast_controller = FilteredIntegrator(gains_fast[j], 0.999, ar1_filter(cutoff_freqs[k], f_loop, "high"), 1/f_loop)
                    sim = Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, 500.0);
                    X_errors[i,j,k] = notched_error_X(sim)
                end
            end
        end

        idx_xmin = argmin(X_errors)
        slow_controller = FilteredIntegrator(gains_slow[idx_xmin[1]], 0.999, ZPKFilter(0, 0, 1), R/f_loop)
        fast_controller = FilteredIntegrator(gains_fast[idx_xmin[2]], 0.999, ar1_filter(cutoff_freqs[idx_xmin[3]], f_loop, "high"), 1/f_loop)
        sim = Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, 500.0)
        search_gain!(sim, "slow")
        search_gain!(sim, "fast")

        idx_xmin_unfiltered = argmin(X_errors_nohpf)
        slow_controller = FilteredIntegrator(gains_slow[idx_xmin_unfiltered[1]], 0.999, ZPKFilter(0, 0, 1), R/f_loop)
        fast_controller = FilteredIntegrator(gains_fast[idx_xmin_unfiltered[2]], 0.999, no_filter, 1/f_loop)
        sim_unfiltered = Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, 500.0)
        search_gain!(sim_unfiltered, "slow")
        search_gain!(sim_unfiltered, "fast")

        push!(opt_xerrs, notched_error_X(sim))
        push!(opt_xerrs_nofilter, notched_error_X(sim_unfiltered))
        push!(opt_xparams, [gains_slow[idx_xmin[1]], gains_fast[idx_xmin[2]], cutoff_freqs[idx_xmin[3]]])
    end
    r0_ncp_labels = [0.6, 1.0, 2.0, 3.0, 6.0]
    plot(r0_ncp_vals, opt_xerrs, xscale=:log10, xticks=(r0_ncp_labels, r0_ncp_labels), label="HPF", xlims=(0.5, 7.0), ylims=(0.5, 1.0), xlabel="NCP r₀ (m)", ylabel="CL error at X (rad)")
    plot!(r0_ncp_vals, opt_xerrs_nofilter, label="No HPF")
end

