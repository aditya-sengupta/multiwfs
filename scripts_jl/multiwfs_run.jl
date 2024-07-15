using multiwfs

loop_cutoffs = 0.0:0.1:40.0

function margins_arr(sys::AOSystem, loop_cutoffs)
    gain_margins, phase_margins = [], []
    for l in loop_cutoffs
        sys.α = exp(-2π * l / sys.f_loop)
        mv = margin(Hol(sys, tf('s')))
        push!(gain_margins, mv.gm[1,1])
        push!(phase_margins, mv.pm[1,1])
    end
    return (gain_margins=gain_margins, phase_margins=phase_margins)
end

base_margins = margins_arr(unfiltered, loop_cutoffs)
low_margins = margins_arr(low_sys, loop_cutoffs)
high_margins = margins_arr(high_sys, loop_cutoffs)

begin
    plot(loop_cutoffs, low_margins.gain_margins, label="LPF gain margin", xlabel="Filter cutoff frequency")
    plot!(loop_cutoffs, high_margins.gain_margins, label="HPF gain margin")
end

begin
    # plot(loop_cutoffs, base_margins.phase_margins)
    plot(loop_cutoffs, low_margins.phase_margins, label="LPF phase margin", xlabel="Filter cutoff frequency")
    plot!(loop_cutoffs, high_margins.phase_margins, label="HPF phase margin")
end

begin
    l = 50.0
    f = (0.1, 500.0)
    high_sys.α = exp(-2π * l / high_sys.f_loop)
    linfreq = range(minimum(f), maximum(f), length=1001)
    linfreq = vcat(-reverse(linfreq), linfreq)
    Hol2plotp = Hol.(Ref(high_sys), linfreq)
    p = plot(real(Hol2plotp), imag(Hol2plotp), xlim=(-1.1,1.1), ylim=(-1.1,1.1), aspect_ratio=:equal, legend=:topright, label="Nyquist plot")
    phasegrid = range(-π, π, length=500)
    xunit, yunit = cos.(phasegrid), sin.(phasegrid)
    plot!(xunit, yunit, ls=:dash, label=nothing)
    vline!([-1/2.5], ls=:dash, label="Gain margin cutoff")
    gain_margin_point = @. abs(imag(Hol2plotp))
    gain_margin_ind = argmin(gain_margin_point)
    scatter!([real(Hol2plotp[gain_margin_ind])], [imag(Hol2plotp[gain_margin_ind])], label="Gain margin point", color=3)
    phase_margin_point = @. abs(real(Hol2plotp)^2 + imag(Hol2plotp)^2 - 1)
    phase_margin_ind = argmin(phase_margin_point)
    plot!([-2,0,-2], [-2,0,2], ls=:dash, label="Phase margin cutoff", color=4)
    scatter!([real(Hol2plotp[phase_margin_ind])], [imag(Hol2plotp[phase_margin_ind])], label="Phase margin point", color=4)
    p
end

