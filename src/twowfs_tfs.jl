using QuadGK

function plant(sT, Cfast, Cslow, R)
    wfs_or_zoh = (1 - exp(-sT)) / sT
    computational_delay = exp(-sT)
    fast_term = Cfast(sT) * computational_delay * wfs_or_zoh^2
    slow_term = Cslow(sT) * computational_delay * wfs_or_zoh * (1 - exp(-sT * R)) / (sT * R)
    return fast_term + slow_term
end

function phi_to_X(z, Cfast, Cslow, R)
    return 1 / (1 + plant(z, Cfast, Cslow, R))
end

function phi_to_Y(z, Cfast, Cslow, R)
    return phi_to_X(z, Cfast, Cslow, R)
end

function Nfast_to_X(z, Cfast, Cslow, R)
    return -(1/z) * Cfast(z) / (1 + plant(z, Cfast, Cslow, R))
end

function Nslow_to_X(z, Cfast, Cslow, R)
    return -(1/z) * Cslow(z) / (1 + plant(z, Cfast, Cslow, R))
end

function Lfast_to_X(z, Cfast, Cslow, R)
    return -(1/z^2) * Cfast(z) / (1 + plant(z, Cfast, Cslow, R))
end

function Lfast_to_Y(z, Cfast, Cslow, R)
    zinv = 1/z
    return (1 + zinv * Cslow(z) * (1/R) * sum(zinv^k for k in 1:R)) / (1 + plant(z, Cfast, Cslow, R))
end

function Lslow_to_X(z, Cfast, Cslow, R)
    return (1 + (1/z^2) * Cfast(z)) / (1 + plant(z, Cfast, Cslow, R))
end

function Lslow_to_Y(z, Cfast, Cslow, R)
    zinv = 1/z
    return -zinv * Cslow(z) * (1/R) * sum(zinv^k for k in 1:R) / (1 + plant(z, Cfast, Cslow, R))
end

function notched_error_X(Cfast, Cslow, R, vk_atm::VonKarman, vk_ncp::VonKarman; f_min=0.1, f_loop=1000.0)
    f_max = f_loop / 2
    f_to_z = f -> exp(2π * im * f / f_loop)
    atm_error_squared = quadgk(f -> psd_von_karman(f, vk_atm) * abs2(phi_to_X(f_to_z(f), Cfast, Cslow, R)), f_min, f_max)[1]
    ncp_error_squared = quadgk(f -> psd_von_karman(f, vk_ncp) * (abs2(Lfast_to_X(f_to_z(f), Cfast, Cslow, R)) + abs2(Lslow_to_X(f_to_z(f), Cfast, Cslow, R))), f_min, f_max)[1]
    return sqrt(atm_error_squared + ncp_error_squared)
end

function notched_error_Y(Cfast, Cslow, R, vk_atm::VonKarman, vk_ncp::VonKarman; f_min=0.1, f_loop=1000.0)
    f_max = f_loop / 2
    f_to_z = f -> exp(2π * im * f / f_loop)
    atm_error_squared = quadgk(f -> psd_von_karman(f, vk_atm) * abs2(phi_to_X(f_to_z(f), Cfast, Cslow, R)), f_min, f_max)[1]
    ncp_error_squared = quadgk(f -> psd_von_karman(f, vk_ncp) * (abs2(Lfast_to_Y(f_to_z(f), Cfast, Cslow, R)) + abs2(Lslow_to_Y(f_to_z(f), Cfast, Cslow, R))), f_min, f_max)[1]
    return sqrt(atm_error_squared + ncp_error_squared)
end

export plant, phi_to_X, phi_to_Y, Nfast_to_X, Nslow_to_X, Lfast_to_X, Lslow_to_X, Lslow_to_Y, Lfast_to_Y, notched_error_X, notched_error_Y