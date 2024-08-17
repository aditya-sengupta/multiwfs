using multiwfs
using Plots
using SciPy
using Optim: optimize, minimizer, Options

function Hzoh_slow(s, slowness=10)
    return (1.0 - exp(-slowness * s)) / (slowness * s)
end

f_loop = 200.0
no_filter = ZPKFilter(0, 0, 1)
sys_test = AOSystem(f_loop, 1.0, 0.3, 0.999, 10, no_filter)

fr = (1e-3:1e-3:0.5) * sys_test.f_loop
sr = 2π * im * fr / sys_test.f_loop
zoh_curve = abs.(Hzoh_slow.(sr))

function to_opt(p)
    approximating_filter = nothing
    try
        approximating_filter = ZPKFilter(signal.cheby1(2, p[1], [p[2]], "lowpass", output="zpk")...)
    catch
        return Inf
    end
    cheb_curve = abs.(transfer_function.(Ref(approximating_filter), sr))
    return sum(abs2, cheb_curve - zoh_curve)
end

res = optimize(to_opt, [0.15, 0.05], SimulatedAnnealing(), Options(iterations=10_000))
minres = minimizer(res)
approximating_filter = ZPKFilter(signal.cheby1(2, minres[1], [minres[2]], "lowpass", output="zpk")...)

begin
    fr = (1e-3:1e-3:0.5) * sys_test.f_loop
    sr = 2π * im * fr / sys_test.f_loop
    plot(fr, abs.(Hzoh_slow.(sr)), xlabel="Frequency (Hz)", label="|Hzoh(f)|")
    plot!(fr, abs.(transfer_function.(Ref(approximating_filter), sr)), label="Approximation")
end