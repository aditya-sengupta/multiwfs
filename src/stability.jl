using NPZ
using Roots
using Optim: optimize

function zero_crossings(f, x, y)
    crossing_regions = findall(x -> x < 0, y[1:end-1] .* y[2:end])
    zero_vals = zeros(length(crossing_regions))
    for (i, crossing_idx) in enumerate(crossing_regions)
        zero_vals[i] = fzero(f, (x[crossing_idx], x[crossing_idx+1]))
    end
    zero_vals
end

"""
Takes in a complex number z and returns its angle relative to -1 in either direction, in degrees.
-1-1im and -1+1im should both return 45.
"""
function angle_relative_to_minus1(z)
    return 180 - abs(rad2deg(angle(z)))
end

function nyquist_and_margins(sim)
    gm, gm_point, pm, pm_point = Inf, nothing, 180, nothing
    oneside_freq = range(0.001, sim.f_loop / 2 + 0.001, length=2001)
    twoside_fr = vcat(-reverse(oneside_freq), oneside_freq)
    twoside_sT = 2π * im * twoside_fr / sim.f_loop
    nyquist_contour = plant.(twoside_sT, Ref(sim))
    imag_axis_crossings = zero_crossings(freq -> imag(plant(2π * im * freq / sim.f_loop, sim)), twoside_fr, imag.(nyquist_contour))
    gm_candidate_points = plant.(imag_axis_crossings * 2π * im / sim.f_loop, Ref(sim))
    if length(gm_candidate_points) > 0
        gm_point = minimum(real(x) for x in gm_candidate_points if (!isnan(x)) & (real(x) > -1))
        gm = -1 / real(gm_point)
    end
    unit_circle_crossings = zero_crossings(freq -> abs2(plant(freq * 2π * im / sim.f_loop, sim)) - 1, twoside_fr, abs2.(nyquist_contour) .- 1)
    pm_candidate_points = plant.(unit_circle_crossings * 2π * im / sim.f_loop, Ref(sim))
    if length(pm_candidate_points) > 0
        pm_point = pm_candidate_points[argmin(angle_relative_to_minus1.(pm_candidate_points))]
        pm = angle_relative_to_minus1(pm_point)
    end
    return nyquist_contour, gm, gm_point, pm, pm_point
end

function margins(sim)
    _, gm, _, pm, _ = nyquist_and_margins(sim)
    return (gm=gm, pm=pm)
end

function is_stable(sim)
    try
        _, gm, _, pm, _ = nyquist_and_margins(sim)
        return is_stable(gm, pm)
    catch
        return false
    end
end

function is_stable(gm, pm)
    return gm > 2.5 && pm >= 45.0
end

function search_gain!(sim, controller_name)
    con = nothing
    if controller_name == "fast"
        con = sim.fast_controller
        con.gain = 1.0
    elseif controller_name == "slow"
        con = sim.slow_controller
        con.gain = 2.0
    end
    gain_min, gain_max = 1e-15, con.gain
    while gain_max - gain_min > 1e-15
        if is_stable(sim)
            gain_min = con.gain
        else
            gain_max = con.gain
        end
        con.gain = (gain_min + gain_max) / 2
    end
    if !is_stable(sim)
        con.gain = con.gain - 1e-15
    end
    con.gain
end

function zero_db_bandwidth(sim)
    abs_etf = f -> abs(phi_to_X(2π * im * f / sim.f_loop, sim))
    upper_limit = optimize(f -> -abs_etf(f), 0.1, 500.0).minimizer
    return find_zero(f -> abs_etf(f) - 1.0, (0.1, upper_limit))
end    

export search_gain!, zero_db_bandwidth, get_nyquist_contour, margins, is_stable