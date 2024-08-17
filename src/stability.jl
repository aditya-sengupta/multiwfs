using NPZ

function nyquist_and_margins(sys; d=0.001)
    f = (0.001, sys.f_loop / 2 + 0.001)
    linfreq = range(minimum(f), maximum(f), length=2001)
    linfreq = vcat(-reverse(linfreq), linfreq)
    nyquist_contour = Hol.(Ref(sys), linfreq)
    gain_margin_points = findall(abs.(imag(nyquist_contour)) .< d)
    @assert length(gain_margin_points) > 0
    gain_margin_ind = argmin(real(nyquist_contour)[gain_margin_points])
    gm = -1 / real(nyquist_contour[gain_margin_points][gain_margin_ind])
    phase_margin_circle_distance = @. abs(abs(nyquist_contour) - 1)
    phase_margin_points = []
    phase_margin_points = findall(phase_margin_circle_distance .< d)
    phase_margin_minus1_distance = @. abs((real(nyquist_contour) + 1)^2 + imag(nyquist_contour)^2)
    phase_margin_ind = argmin(phase_margin_minus1_distance[phase_margin_points])
    pm = abs(rad2deg((angle(nyquist_contour[phase_margin_points][phase_margin_ind]) - π) % (2π)))
    return nyquist_contour, gm, gain_margin_points[gain_margin_ind], pm, phase_margin_points[phase_margin_ind]
end

function is_stable(sys)
    try
        _, gm, _, pm, _ = nyquist_and_margins(sys)
        return is_stable(gm, pm)
    catch
        return false
    end
end

function is_stable(gm, pm)
    return gm > 2.5 && pm >= 45.0
end

function search_gain!(sys)
    sys.gain = 1.0
    gain_min, gain_max = 1e-15, 1.0
    while gain_max - gain_min > 1e-15
        if is_stable(sys)
            gain_min = sys.gain
        else
            gain_max = sys.gain
        end
        sys.gain = (gain_min + gain_max) / 2
    end
    if !is_stable(sys)
        sys.gain = sys.gain - 1e-15
    end
    sys.gain
end

function zero_db_bandwidth(sys)
    try
        # try to solve via root-finding
        return find_zero(f -> abs(Hol(sys, f)) - 1.0, (0.1, 500.0))
    catch
        # fall back to grid evaluation
        f = 0.1:0.1:500.0
        abs_Hol_val = abs.(Hol.(Ref(sys), f))
        fstart = argmax(abs_Hol_val)
        fend = findlast(abs_Hol_val .<= 1.0)
        return f[fstart:fend][findfirst(abs_Hol_val[fstart:fend] .<= 1.0)]
    end
end    

function ar1_gain_map(sys, filter_type; f_cutoffs = 0.1:0.1:100.0, delays = 0.0:0.1:1.0, save=true)
    gain_map = zeros(length(f_cutoffs), length(delays));
    @showprogress @threads for (i, fc) in collect(enumerate(f_cutoffs))
        for (j, d) in enumerate(delays)
            tsys = AOSystem(sys.f_loop, d, sys.gain, sys.leak, sys.fpf, ar1_filter(fc, sys.f_loop, filter_type))
            search_gain!(tsys)
            gain_map[i,j] = tsys.gain
        end
    end
    if save
        npzwrite("data/gainmap_loopfreq_$(sys.f_loop)_ftype_$filter_type.npy", gain_map)
    end
    gain_map
end

export ar1_gain_map, search_gain!, zero_db_bandwidth