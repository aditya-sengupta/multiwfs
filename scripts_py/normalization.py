# %%
import numpy as np
from scipy.signal import welch

time_series = np.load("../data/olt/olt_atm.npy")[:1000]
variance = np.var(time_series)

# Compute the power spectral density using Welch's method
frequencies, power_spectrum = welch(time_series, fs=1.0, scaling='density')

# Compute the area under the power spectrum
area_under_power_spectrum = np.sum(power_spectrum) * (frequencies[1] - frequencies[0])

print(f"Variance of the time-series: {variance}")
print(f"Area under the power spectrum: {area_under_power_spectrum}")
# %%
