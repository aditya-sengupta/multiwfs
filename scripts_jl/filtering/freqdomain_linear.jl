include("filter_setup.jl")

begin
    f_cutoff = 15.0
    # parameters are:
    # f_loop, frame_delay, gain, leak, fpf, filter
    f_loop = 1000.0
    sys_slow = AOSystem(f_loop, 1.0, 1.4, 0.999, 1, no_filter)
    sys_fast = AOSystem(f_loop, 0.1, 0.4, 0.999, 10, ar1_high)
    #search_gain!(sys_low)
    #search_gain!(sys_high)

    f = exp10.(-2:0.01:log10(f_loop/2))

    Hol_vals = Hol.(Ref(sys_slow), f) .+ Hol.(Ref(sys_fast), f)
    Hrej_vals = @. 1 / (1 + Hol_vals)

    plot(f, abs2.(Hrej_vals), xscale=:log10, yscale=:log10, xticks=[1e-2, 1e-1, 1e0, 1e1, 1e2], xlabel="Frequency (Hz)", ylabel="Error transfer function", label=nothing)
end

