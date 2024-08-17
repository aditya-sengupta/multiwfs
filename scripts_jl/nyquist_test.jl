using multiwfs
include("filter_setup.jl")

f_loop = 200.0
f_cutoff = 3.0

sys_high = AOSystem(
    f_loop, 1.0, 0.3, 0.999, 10,
    cheb2_high
)
search_gain!(sys_high)
nyquist_plot(sys_high; d=0.001)

