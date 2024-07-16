using multiwfs
using Plots
using Roots

sys = AOSystem(1000.0, 1.0, 0.01, 0.999, 10, "high", 80.0)
search_gain!(sys)
nyquist_plot(sys)
zero_db_bandwidth(sys)