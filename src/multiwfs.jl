module multiwfs
    using Plots
    using StatsBase: median
    using ProgressMeter
    using Base.Threads
    using NPZ

    function f2s(f)
        return 1im * 2.0 * π * f
    end

    mutable struct AOSystem
        f_loop::Float64
        Ts::Float64
        frame_delay::Float64
        τ::Float64
        gain::Float64
        leak::Float64
        fpf::Int64
        filter_type::String
        filter_cutoff::Float64
        α::Float64

        function AOSystem(f_loop, frame_delay, gain, leak, fpf, filter_type, filter_cutoff)
            Ts = 1 / f_loop
            frame_delay = floor(frame_delay) + round((frame_delay-floor(frame_delay))*fpf)/fpf
            τ = frame_delay * Ts
            α = exp(-2π * filter_cutoff / f_loop)
            new(f_loop, Ts, frame_delay, τ, gain, leak, fpf, filter_type, filter_cutoff, α)
        end
    end

    function update_filter_cutoff!(sys, filter_cutoff)
        sys.filter_cutoff = filter_cutoff
        sys.α = exp(-2π * filter_cutoff / sys.f_loop)
    end

    function update_frame_delay!(sys, frame_delay)
        sys.frame_delay = floor(frame_delay) + round((frame_delay-floor(frame_delay))*sys.fpf)/sys.fpf
        sys.τ = sys.frame_delay * sys.Ts
    end

    function Hwfs(ao::AOSystem, s)
        return (1.0 - exp(-ao.Ts * s)) / (ao.Ts * s)
    end

    function Hzoh(ao::AOSystem, s)
        return Hwfs(ao, s)
    end

    function Hlag(ao::AOSystem, s)
        return exp(-ao.τ * s)
    end

    function Hint(ao::AOSystem, s)
        return ao.gain / (1.0 - exp(-ao.Ts * s)) # integrator
    end

    function Hlint(ao::AOSystem, s)
        return ao.gain / (1.0 - ao.leak * exp(-ao.Ts * s)) # leaky integrator
    end

    function Hcont(ao::AOSystem, s)
        return Hlint(ao, s)
    end

    function Hfilter(ao::AOSystem, s)
        if ao.filter_type == "high"
            return ao.α * (1 - exp(-s / ao.f_loop)) / (1 - ao.α * exp(-s / ao.f_loop))
        elseif ao.filter_type == "low"
            return (1 - ao.α) / (1 - ao.α * exp(-s / ao.f_loop))
        else
            return 1
        end
    end

    function Hol(ao::AOSystem, s::Complex)
        return Hwfs(ao, s) * Hlag(ao, s) * Hcont(ao, s) * Hfilter(ao, s) * Hzoh(ao, s)
    end

    function Hol_unfiltered(ao::AOSystem, f::Real)
        s = f2s(f)
        return Hwfs(ao, s) * Hlag(ao, s) * Hcont(ao, s) * Hzoh(ao, s)
    end

    function Hol(ao::AOSystem, f)
        s = f2s(f)
        return Hwfs(ao, s) * Hlag(ao, s) * Hcont(ao, s) * Hfilter(ao, s) * Hzoh(ao, s)
    end

    function Hrej(ao::AOSystem, f)
        return 1.0 / (1.0 + Hol(ao, f))
    end

    function Hcl(ao::AOSystem, f)
        Hol_val = Hol(ao, f)
        return Hol_val / (1.0 + Hol_val)
    end

    function Hn(ao::AOSystem, f)
        return Hcl(ao, f) / Hwfs(ao, f2s(f))
    end

    function nyquist_and_margins(sys)
        f = (0.001, sys.f_loop / 2 + 0.001)
        linfreq = range(minimum(f), maximum(f), length=2001)
        linfreq = vcat(-reverse(linfreq), linfreq)
        nyquist_contour = Hol.(Ref(sys), linfreq)
        gain_margin_points = findall(abs.(imag(nyquist_contour)) .< 0.001)
        @assert length(gain_margin_points) > 0
        gain_margin_ind = argmin(real(nyquist_contour)[gain_margin_points])
        gm = -1 / real(nyquist_contour[gain_margin_points][gain_margin_ind])
        phase_margin_points = @. abs(real(nyquist_contour)^2 + imag(nyquist_contour)^2 - 1)
        phase_margin_ind = argmin(phase_margin_points)
        pm = abs(rad2deg((angle(nyquist_contour[phase_margin_ind]) - π) % (2π)))
        return nyquist_contour, gm, gain_margin_points[gain_margin_ind], pm, phase_margin_ind
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

    function nyquist_plot(sys)
        success = :green
        nyquist_contour, gm, gain_margin_ind, pm, phase_margin_ind = nyquist_and_margins(sys)
        p = plot(real(nyquist_contour), imag(nyquist_contour), xlim=(-1.1,1.1), ylim=(-1.1,1.1), aspect_ratio=:equal, legend=:outertopright, label="Nyquist plot")
        phasegrid = range(-π, π, length=500)
        xunit, yunit = cos.(phasegrid), sin.(phasegrid)
        vline!([-1/2.5], ls=:dash, label="Gain margin cutoff", color=:grey)
        if !is_stable(gm, pm)
            success = :red
        end
        scatter!([real(nyquist_contour[gain_margin_ind])], [imag(nyquist_contour[gain_margin_ind])], label="Gain margin = $(round(gm, digits=2))", color=:grey)
        plot!([-2,0,-2], [-2,0,2], ls=:dash, label="Phase margin cutoff", color=4)
        plot!(xunit, yunit, ls=:dash, label=nothing, color=success)
        scatter!([real(nyquist_contour[phase_margin_ind])], [imag(nyquist_contour[phase_margin_ind])], label="Phase margin = $(round(pm, digits=2))", color=4, title="$(uppercase((sys.filter_type[1])))PF, f_cutoff=$(sys.filter_cutoff) Hz, delay=$(sys.frame_delay)")
        p
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

    function gain_map(sys; f_cutoffs = 0.1:0.1:100.0, delays = 0.0:0.1:1.0, save=true)
        gain_map = zeros(length(f_cutoffs), length(delays));
        @showprogress @threads for (i, fc) in collect(enumerate(f_cutoffs))
            for (j, d) in enumerate(delays)
                tsys = AOSystem(sys.f_loop, d, sys.gain, sys.leak, sys.fpf, sys.filter_type, fc)
                search_gain!(tsys)
                gain_map[i,j] = tsys.gain
            end
        end
        if save
            npzwrite("data/gainmap_loopfreq_$(sys.f_loop)_ftype_$(sys.filter_type).npy", gain_map)
        end
        gain_map
    end    

    function psd(x)
        noverlap = 2^7# 2^(Int(floor(log2(length(x)/4)))-2)
        n = div(length(x), 8)
        return welch_pgram(x, n, noverlap; fs=f_loop)
    end
    
    function plot_psd_p!(f, p; kwargs...)
        plot!(f, p, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="Power"; kwargs...)
    end
    
    function integrator_control(open_loop, gain, leak, update_every; hpf_gain=0.0, delay_frames=1)
        closed_loop = zeros(length(open_loop))
        closed_loop[1] = open_loop[1]
        command = 0.0
        average_buffer = []
        hpf_val = 0.0
        for i in 2:N
            push!(average_buffer, closed_loop[i-1])
            if i > 2 + delay_frames
                hpf_val = sys_high.α * hpf_val + sys_high.α * (closed_loop[i-1-delay_frames] - closed_loop[i-2-delay_frames])
            end
            if i % update_every == 0
                command = leak * command - gain * mean(average_buffer)
                average_buffer = []
            end
            # hpf is stuck to update_every = 1
            command -= hpf_gain * hpf_val 
            closed_loop[i] = open_loop[i] + command
        end
        closed_loop
    end

    export AOSystem, nyquist_and_margins, nyquist_plot, is_stable, search_gain!, zero_db_bandwidth, update_filter_cutoff!, update_frame_delay!, gain_map, psd, plot_psd_p!, integrator_control, Hrej, Hol_unfiltered, Hol

end # module multiwfs
