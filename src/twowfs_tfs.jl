using QuadGK
using Base.Meta: parse

function plant(sT, Cfast, Cslow, R, dm_slow_rate=true)
    wfs_or_zoh = (1 - exp(-sT)) / sT
    computational_delay = exp(-sT)
    fast_term = Cfast(sT) * computational_delay * wfs_or_zoh^2
    if dm_slow_rate
        slow_term = Cslow(sT) * computational_delay * ((1 - exp(-sT * R)) / (sT * R))^2
    else
        slow_term = Cslow(sT) * computational_delay * wfs_or_zoh * ((1 - exp(-sT * R)) / (sT * R))

    end
    return fast_term + slow_term
end

function phi_to_X(sT, Cfast, Cslow, R, dm_slow_rate=true)
    return 1 / (1 + plant(sT, Cfast, Cslow, R, dm_slow_rate))
end

function phi_to_Y(sT, Cfast, Cslow, R)
    return phi_to_X(sT, Cfast, Cslow, R)
end

function Lfast_to_X(sT, Cfast, Cslow, R)
    esT = exp(-sT)
    return -((1 - esT)/sT)^2 * esT * Cfast(sT) / (1 + plant(sT, Cfast, Cslow, R))
end

function Lfast_to_Y(sT, Cfast, Cslow, R)
    esT = exp(-sT)
    return (1 + ((1 - esT) / sT) * esT * Cslow(sT) * ((1 - exp(-sT*R)) / (sT * R))) / (1 + plant(sT, Cfast, Cslow, R))
end

function Lslow_to_X(sT, Cfast, Cslow, R)
    esT = exp(-sT)
    return (1 + ((1 - esT) / sT)^2 * esT * Cfast(sT)) / (1 + plant(sT, Cfast, Cslow, R))
end

function Lslow_to_Y(sT, Cfast, Cslow, R)
    esT = exp(-sT)
    return -((1 - esT) / sT) * esT * Cslow(sT) * (1 - exp(-sT*R)) / (sT * R) / (1 + plant(sT, Cfast, Cslow, R))
end

function noise_tf_prefactor(sT, Cfast, Cslow, R)
    esT = exp(-sT)
    half_delay = (1 - esT) / sT
    numerator = -half_delay * esT
    denominator = 1 + plant(sT, Cfast, Cslow, R)
    return numerator / denominator
end

function Nslow_to_X(sT, Cfast, Cslow, R)
    if imag(sT) < π / R
        return noise_tf_prefactor(sT, Cfast, Cslow, R) * Cslow(sT)
    else
        return 0.0
    end
end

function Nslow_to_Y(sT, Cfast, Cslow, R)
    return Nslow_to_X(sT, Cfast, Cslow, R)
end

function Nfast_to_X(sT, Cfast, Cslow, R)
    return noise_tf_prefactor(sT, Cfast, Cslow, R) * Cfast(sT)
end

function Nfast_to_Y(sT, Cfast, Cslow, R)
    return Nfast_to_X(sT, Cfast, Cslow, R)
end

function atm_error_at_f_X(f, Cfast, Cslow, R, vk_atm::VonKarman, vk_ncp::VonKarman, noise_normalization, f_loop)
    sT = 2π * im * f / f_loop
    psd_von_karman(f, vk_atm) * abs2(phi_to_X(sT, Cfast, Cslow, R))
end

function ncp_error_at_f_X(f, Cfast, Cslow, R, vk_atm::VonKarman, vk_ncp::VonKarman, noise_normalization, f_loop)
    sT = 2π * im * f / f_loop
    psd_von_karman(f, vk_ncp) * (abs2(Lfast_to_X(sT, Cfast, Cslow, R)) + abs2(Lslow_to_X(sT, Cfast, Cslow, R)))
end

function noise_error_at_f_X(f, Cfast, Cslow, R, vk_atm::VonKarman, vk_ncp::VonKarman, noise_normalization, f_loop)
    sT = 2π * im * f / f_loop
    noise_normalization * (abs2(Nfast_to_X(sT, Cfast, Cslow, R)) + abs2(Nslow_to_X(sT, Cfast, Cslow, R)))
end

function error_at_f_X(f, Cfast, Cslow, R, vk_atm::VonKarman, vk_ncp::VonKarman, noise_normalization, f_loop)
    return atm_error_at_f_X(f, Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization, f_loop) + ncp_error_at_f_X(f, Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization, f_loop) + noise_error_at_f_X(f, Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization, f_loop)
end

function atm_error_at_f_Y(f, Cfast, Cslow, R, vk_atm::VonKarman, vk_ncp::VonKarman, noise_normalization, f_loop)
    sT = 2π * im * f / f_loop
    psd_von_karman(f, vk_atm) * abs2(phi_to_Y(sT, Cfast, Cslow, R))
end

function ncp_error_at_f_Y(f, Cfast, Cslow, R, vk_atm::VonKarman, vk_ncp::VonKarman, noise_normalization, f_loop)
    sT = 2π * im * f / f_loop
    psd_von_karman(f, vk_ncp) * (abs2(Lfast_to_Y(sT, Cfast, Cslow, R)) + abs2(Lslow_to_Y(sT, Cfast, Cslow, R)))
end

function noise_error_at_f_Y(f, Cfast, Cslow, R, vk_atm::VonKarman, vk_ncp::VonKarman, noise_normalization, f_loop)
    sT = 2π * im * f / f_loop
    noise_normalization * (abs2(Nfast_to_Y(sT, Cfast, Cslow, R)) + abs2(Nslow_to_Y(sT, Cfast, Cslow, R)))
end

function error_at_f_Y(f, Cfast, Cslow, R, vk_atm::VonKarman, vk_ncp::VonKarman, noise_normalization, f_loop)
    return atm_error_at_f_Y(f, Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization, f_loop) + ncp_error_at_f_Y(f, Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization, f_loop) + noise_error_at_f_Y(f, Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization, f_loop)
end

function notched_error_X(Cfast, Cslow, R, vk_atm::VonKarman, vk_ncp::VonKarman, f_noise_crossover; f_min=0.1, f_loop=1000.0)
    f_max = f_loop / 2
    noise_normalization = psd_von_karman(f_noise_crossover, vk_atm)
    return sqrt(
        quadgk(f -> error_at_f_X(f, Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization, f_loop), f_min, f_max)[1]
    )
end

function notched_error_Y(Cfast, Cslow, R, vk_atm::VonKarman, vk_ncp::VonKarman, f_noise_crossover; f_min=0.1, f_loop=1000.0)
    f_max = f_loop / 2
    noise_normalization = psd_von_karman(f_noise_crossover, vk_atm)
    return sqrt(
        quadgk(f -> error_at_f_Y(f, Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization, f_loop), f_min, f_max)[1]
    )
end

export plant, phi_to_X, phi_to_Y, Nfast_to_X, Nslow_to_X, Nfast_to_Y, Nslow_to_Y, Lfast_to_X, Lslow_to_X, Lslow_to_Y, Lfast_to_Y, error_at_f_X, atm_error_at_f_X, ncp_error_at_f_X, noise_error_at_f_X, atm_error_at_f_Y, ncp_error_at_f_Y, noise_error_at_f_Y, error_at_f_Y, notched_error_X, notched_error_Y