using multiwfs
using multiwfs: Hol
using Plots
using Roots

sys = AOSystem(1000.0, 1.0, 0.01, 0.999, 10, "high", 80.0)
search_gain!(sys)
nyquist_plot(sys)
plot(f -> abs(Hol(sys, f)), 0.1:0.1:500.0)
zero_db_bandwidth(sys)