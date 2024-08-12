using multiwfs
using Plots

function circle_contour(h, k, r, θmin=0, θmax=2π)
    θ = LinRange(θmin, θmax, 500)
    h .+ r .* sin.(θ), k .+ r .* cos.(θ)
end

sys = AOSystem(1000.0, 1.0, 0.01, 0.999, 10, "low", 80.0)
# gain_map(sys)
search_gain!(sys)
sys.gain = 0.45
nyquist_plot(sys)

nyquist_contour, gm, gain_margin_ind, pm, phase_margin_ind = nyquist_and_margins(sys)

etf_contour = @. 1 / (1 + nyquist_contour)

begin
    plot(real.(etf_contour), imag.(etf_contour), aspect_ratio=:equal, ylim=(-2,2), xlim=(-2,2))
    #plot!(circle_contour(1/2, 1/2, sqrt(2))..., ls=:dash)
    #plot!(circle_contour(1/2, -1/2, sqrt(2))..., ls=:dash)
    plot!(circle_contour(5/3, 0, 5/3)..., ls=:dash)
    plot!(circle_contour(0, 0, 1.31)..., ls=:dash)
end