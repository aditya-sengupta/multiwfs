using multiwfs
using JLD2
using Plots
using Colors
using FITSIO
using DataInterpolations
using QuadGK

fast_freq, fast_atmplusnoise, fast_ncpplusnoise, fast_noise = eachcol(read(FITS("data/data_fast_psds.fits")[1]))
slow_freq, slow_atmplusnoise, slow_ncpplusnoise, slow_noise = eachcol(read(FITS("data/data_slow_psds.fits")[1]))
optpars = load("data/all_opt.jld2")
threem_pars = optpars["(1.0, 200.0)"]
vk_ncp = VonKarman(0.001, 0.25 * 1.0^(-5/3))
f_crossover = 200.0
sim_one = simgen_ichpf(threem_pars[1][1], 0.0, 0.0, vk_ncp, f_crossover; leak=0.999)
fast_atmplusnoise_interp = CubicSpline(fast_atmplusnoise, fast_freq)

begin
    plot(fast_freq, fast_atmplusnoise, xscale=:log10, yscale=:log10, color=1, label="atm, time", yticks=[1e-6, 1e-4, 1e-2, 1e0], xlabel="Frequency (Hz)", ylabel="Power spectral density (rad²/Hz)")
    plot!(fast_freq, psd_von_karman.(fast_freq, Ref(sim_one.vk_atm)), color=1, ls=:dash, label="atm, freq")
    plot!(fast_freq, fast_ncpplusnoise, xscale=:log10, yscale=:log10, color=2, label="ncp, time")
    plot!(fast_freq,psd_von_karman.(fast_freq, Ref(sim_one.vk_ncp)), color=2, ls=:dash, label="ncp, freq")
    plot!(fast_freq, fast_noise, xscale=:log10, yscale=:log10, color=3, label="noise, time")
    plot!(fast_freq, repeat([psd_von_karman(f_crossover, sim_one.vk_atm)], length(fast_freq)), color=3, ls=:dash, label="noise, freq")
end

function error_at_f_X_olpsds(f, sim::Simulation, ncp=true)
    fast_ncpplusnoise_interp = CubicSpline(psd_von_karman.(fast_freq, Ref(sim.vk_ncp)), fast_freq)
    slow_ncpplusnoise_interp = CubicSpline(psd_von_karman.(slow_freq, Ref(sim.vk_ncp)), slow_freq)
    sT = 2π * im * f / sim.f_loop
    err = 0.0
    err += fast_atmplusnoise_interp(f) * abs2(phi_to_X(sT, sim))
    if ncp
        err += fast_ncpplusnoise_interp(f) * abs2(Lfast_to_X(sT, sim))
        if f <= last(slow_freq)
            err += slow_ncpplusnoise_interp(f) * abs2(Lslow_to_X(sT, sim))
        end
    end
    err += sim.noise_normalization * abs2(Nfast_to_X(sT, sim))
    if f <= last(slow_freq)
        err += sim.noise_normalization * abs2(Nslow_to_X(sT, sim))
    end
    err
end

function notched_error_X_olpsds(sim, ncp=true)
    if is_stable(sim)
        return sqrt(quadgk(f -> error_at_f_X_olpsds(f, sim, ncp), sim.f_min_cost, last(fast_freq))[1])
    else
        return Inf
    end
end

notched_error_X(sim_one)
notched_error_X_olpsds(sim_one)

notched_error_X_noncp(sim_one)
notched_error_X_olpsds(sim_one, false)

begin
    xerrs = Dict()
    r0_ncps = vcat(0.6:0.2:2.0, 3:6)
    for r0_ncp in r0_ncps
        vk_ncp = VonKarman(0.001, 0.25 * (r0_ncp)^(-5/3))
        sim_one = simgen_ichpf(optpars["($r0_ncp, $f_crossover)"][1][1], 0.0, 0.0, vk_ncp, f_crossover; leak=0.999)
        sim_two = simgen_ichpf(optpars["($r0_ncp, $f_crossover)"][3]..., vk_ncp, f_crossover; leak=0.999)
        xerrs["($r0_ncp, $f_crossover) one"] = notched_error_X_olpsds(sim_one)
        xerrs["($r0_ncp, $f_crossover) one noncp"] = notched_error_X_olpsds(sim_one, false)
        xerrs["($r0_ncp, $f_crossover) two"] = notched_error_X_olpsds(sim_two)
        xerrs["($r0_ncp, $f_crossover) two noncp"] = notched_error_X_olpsds(sim_two, false)
        #xerrs["($r0_ncp, $f_crossover) one"] = notched_error_X(sim_one)
        #xerrs["($r0_ncp, $f_crossover) one noncp"] = notched_error_X_noncp(sim_one)
        #xerrs["($r0_ncp, $f_crossover) two"] = notched_error_X(sim_two)
        #xerrs["($r0_ncp, $f_crossover) two noncp"] = notched_error_X_noncp(sim_two)
    end
    
    p = Plots.plot(xlabel="r₀ NCP (m)", ylabel="CL error at X (rad)", xscale=:log10, xticks=([0.6, 1.0, 1.4, 2.0, 3.0, 4.0, 6.0], [0.6, 1.0, 1.4, 2.0, 3.0, 4.0, 6.0]), ylims=(0.0, 1.5), suptitle="Hairdryer plot, 200 Hz, AO simulator OL PSDs")
    Plots.plot!(r0_ncps, [xerrs["($r0_ncp, $f_crossover) one noncp"] for r0_ncp in r0_ncps], msw=0, label="Main: Atm only", color=RGBA(1/255,115/255,178/255,255/255), lw=2)
    Plots.plot!(r0_ncps, [xerrs["($r0_ncp, $f_crossover) two noncp"] for r0_ncp in r0_ncps], msw=0, label="Both: Atm only", color=RGBA(222/255,143/255,5/255,255/255), lw=2)
    Plots.plot!(r0_ncps, [xerrs["($r0_ncp, $f_crossover) one"] for r0_ncp in r0_ncps], msw=0, label="Main: Atm + NCP", color=RGBA(2/255,158/255,115/255,255/255), lw=2)
    Plots.plot!(r0_ncps, [xerrs["($r0_ncp, $f_crossover) two"] for r0_ncp in r0_ncps], msw=0, label="Both: Atm + NCP", color=RGBA(213/255,94/255,0/255,255/255), lw=2)
    Plots.plot(p)
end
