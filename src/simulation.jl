using QuadGK
using Base.Meta: parse

struct Simulation{C1,C2}
    f_loop::Float64
    fast_controller::C1
    slow_controller::C2
    fast_plant::Plant
    slow_plant::Plant
    fr::Vector{Float64}
    sr::Vector{ComplexF64}
    sT::Vector{ComplexF64}
    R::Int64
    vk_atm::VonKarman
    vk_ncp::VonKarman
    f_noise_crossover::Float64
    noise_normalization::Float64
    f_min_cost::Float64
    f_max::Float64

    function Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover; f_min_eval=0.01, f_min_cost=1.0, n_freq_points=470, frame_delay=1.0, fpf=1, f_max=nothing)
        if isnothing(f_max)
            f_max = f_loop / 2
        end
        log_f_min_eval = log10(f_min_eval)
        log_f_max = log10(f_loop/2)
        log_f_step = (log_f_max - log_f_min_eval) / (n_freq_points - 1)
        fr = exp10.(log_f_min_eval:log_f_step:log_f_max)
        sr = 2π .* im * fr
        sT = sr .* f_loop
        noise_normalization = psd_von_karman(f_noise_crossover, vk_atm)
        fast_plant = Plant(f_loop, frame_delay, fpf)
        slow_plant = Plant(f_loop / R, frame_delay / R, fpf)

        new{typeof(fast_controller),typeof(slow_controller)}(f_loop, fast_controller, slow_controller, fast_plant, slow_plant, fr, sr, sT, R, vk_atm, vk_ncp, f_noise_crossover, noise_normalization, f_min_cost, f_max)
    end
end

function plant(sT, sim::Simulation, dm_slow_rate=true)
    wfs_or_zoh = (1 - exp(-sT)) / sT
    computational_delay = exp(-sT)
    fast_term = transfer_function(sim.fast_controller, sT) * computational_delay * wfs_or_zoh^2
    if dm_slow_rate
        slow_term = transfer_function(sim.slow_controller, sT) * computational_delay * ((1 - exp(-sT * sim.R)) / (sT * sim.R))^2
    else
        slow_term = transfer_function(sim.slow_controller, sT) * computational_delay * wfs_or_zoh * ((1 - exp(-sT * sim.R)) / (sT * sim.R))

    end
    return fast_term + slow_term
end

function phi_to_X(sT, sim::Simulation, dm_slow_rate=true)
    return 1 / (1 + plant(sT, sim, dm_slow_rate))
end

phi_to_Y = phi_to_X

function Lfast_to_X(sT, sim::Simulation)
    esT = exp(-sT)
    return -((1 - esT)/sT)^2 * esT * transfer_function(sim.fast_controller, sT) / (1 + plant(sT, sim))
end

function Lfast_to_Y(sT, sim::Simulation)
    esT = exp(-sT)
    return (1 + ((1 - esT) / sT) * esT * transfer_function(sim.slow_controller, sT) * ((1 - exp(-sT*sim.R)) / (sT * sim.R))) / (1 + plant(sT, sim))
end

function Lslow_to_X(sT, sim::Simulation)
    esT = exp(-sT)
    return (1 + ((1 - esT) / sT)^2 * esT * transfer_function(sim.fast_controller, sT)) / (1 + plant(sT, sim))
end

function Lslow_to_Y(sT, sim::Simulation)
    esT = exp(-sT)
    return -((1 - esT) / sT) * esT * transfer_function(sim.slow_controller, sT) * (1 - exp(-sT*sim.R)) / (sT * sim.R) / (1 + plant(sT, sim))
end

function noise_tf_prefactor(sT, sim::Simulation)
    esT = exp(-sT)
    half_delay = (1 - esT) / sT
    numerator = -half_delay * esT
    denominator = 1 + plant(sT, sim)
    return numerator / denominator
end

function Nslow_to_X(sT, sim::Simulation)
    if imag(sT) < π / sim.R
        return noise_tf_prefactor(sT, sim) * transfer_function(sim.slow_controller, sT)
    else
        return 0.0
    end
end

Nslow_to_Y = Nslow_to_X

function Nfast_to_X(sT, sim::Simulation)
    return noise_tf_prefactor(sT, sim) * transfer_function(sim.fast_controller, sT)
end

Nfast_to_Y = Nfast_to_X

function atm_error_at_f_X(f, sim::Simulation)
    sT = 2π * im * f / sim.f_loop
    psd_von_karman(f, sim.vk_atm) * abs2(phi_to_X(sT, sim))
end

function ncp_error_at_f_X(f, sim::Simulation)
    sT = 2π * im * f / sim.f_loop
    psd_von_karman(f, sim.vk_ncp) * (abs2(Lfast_to_X(sT, sim)) + abs2(Lslow_to_X(sT, sim)))
end

function noise_error_at_f_X(f, sim::Simulation)
    sT = 2π * im * f / sim.f_loop
    sim.noise_normalization * (abs2(Nfast_to_X(sT, sim)) + abs2(Nslow_to_X(sT, sim)))
end

function error_at_f_X(f, sim::Simulation)
    return atm_error_at_f_X(f, sim) + ncp_error_at_f_X(f, sim) + noise_error_at_f_X(f, sim)
end

function atm_error_at_f_Y(f, sim::Simulation)
    sT = 2π * im * f / sim.f_loop
    psd_von_karman(f, sim.vk_atm) * abs2(phi_to_Y(sT, sim))
end

function ncp_error_at_f_Y(f, sim::Simulation)
    sT = 2π * im * f / sim.f_loop
    psd_von_karman(f, sim.vk_ncp) * (abs2(Lfast_to_Y(sT, sim)) + abs2(Lslow_to_Y(sT, sim)))
end

function noise_error_at_f_Y(f, sim::Simulation)
    sT = 2π * im * f / sim.f_loop
    sim.noise_normalization * (abs2(Nfast_to_Y(sT, sim)) + abs2(Nslow_to_Y(sT, sim)))
end

function error_at_f_Y(f, sim::Simulation)
    return atm_error_at_f_Y(f, sim) + ncp_error_at_f_Y(f, sim) + noise_error_at_f_Y(f, sim)
end

function notched_error_X(sim::Simulation)
    return sqrt(
        quadgk(f -> error_at_f_X(f, sim), sim.f_min_cost, sim.f_max)[1]
    )
end

function notched_error_Y(sim::Simulation)
    return sqrt(
        quadgk(f -> error_at_f_Y(f, sim), sim.f_min_cost, sim.f_max)[1]
    )
end

export Simulation, plant, phi_to_X, phi_to_Y, Nfast_to_X, Nslow_to_X, Nfast_to_Y, Nslow_to_Y, Lfast_to_X, Lslow_to_X, Lslow_to_Y, Lfast_to_Y, error_at_f_X, atm_error_at_f_X, ncp_error_at_f_X, noise_error_at_f_X, atm_error_at_f_Y, ncp_error_at_f_Y, noise_error_at_f_Y, error_at_f_Y, notched_error_X, notched_error_Y