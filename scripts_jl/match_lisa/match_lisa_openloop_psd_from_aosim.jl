using multiwfs
using FITSIO
using Plots

fast_freq, fast_atmplusnoise, fast_ncpplusnoise, fast_noise = eachcol(read(FITS("data/data_fast_psds.fits")[1]))
slow_freq, slow_atmplusnoise, slow_ncpplusnoise, slow_noise = eachcol(read(FITS("data/data_slow_psds.fits")[1]))

begin
    plot(fast_freq, fast_atmplusnoise, xscale=:log10, yscale=:log10, color=:orange, label="Fast WFS atm+noise", lw=2)
    plot!(slow_freq, slow_atmplusnoise, xscale=:log10, yscale=:log10, color=:red, ls=:dash, label="Slow WFS atm+noise", lw=2)

    plot!(fast_freq, fast_ncpplusnoise, xscale=:log10, yscale=:log10, color=:mediumseagreen, label="Fast WFS NCP+noise", lw=2)
    plot!(slow_freq, slow_ncpplusnoise, xscale=:log10, yscale=:log10, color=:darkgreen, ls=:dash, label="Slow WFS NCP+noise", lw=2)

    plot!(fast_freq, fast_noise, xscale=:log10, yscale=:log10, color=:royalblue1, label="Fast WFS noise", lw=2)
    plot!(slow_freq, slow_noise, xscale=:log10, yscale=:log10, color=:darkblue, ls=:dash, label="Slow WFS noise", lw=2)
end