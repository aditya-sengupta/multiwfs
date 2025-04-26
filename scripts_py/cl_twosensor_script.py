# %%
import numpy as np
from matplotlib import pyplot as plt
from cl_twosensor_withlqg import run_cl_2sensor
from lisa_psd import genAvgPerUnb

# %%
atm = np.load("../data/olt/olt_atm.npy")
ncp_fast = np.load("../data/olt/olt_fastncp.npy")
ncp_slow = np.load("../data/olt/olt_slowncp.npy")
noise_fast = np.load("../data/olt/olt_fastnoise.npy")
noise_slow = np.load("../data/olt/olt_slownoise.npy")

# %%
for (ts, n) in zip([atm, ncp_fast, noise_fast], ["atm", "ncp", "noise"]):
    plt.loglog(*genpsd(ts, 1e-3), label=n)
plt.xlabel("Frequency (Hz)")
plt.ylabel(r"Power (rad${}^2$/Hz)")
plt.legend()
# %%
output_ts,  \
meas_slow_ts, \
meas_fast_ts, \
meas_fast_af_ts, \
slowphase_ts= \
    run_cl_2sensor(atm*0, 10, delayfr=1, \
                    slow_gain=1.4, fast_gain=0.4, integ=0.995,
                    ar1hp_coeff=0.91,
                    ncptimeseries=ncp_fast*0,
                    slowncptimeseries=ncp_slow*0,
                    mnoisetimeseries=noise_fast*0,
                    slowmnoisetimeseries=noise_slow
                    )
    

psd_noise_slow = genAvgPerUnb(noise_slow, 2048)[1:1025]
psd_cl = genAvgPerUnb(slowphase_ts, 2048)[1:1025]
freq = np.linspace(0.0, 500.0, 1025)[1:]

#plt.loglog(freq, psd_cl / psd_atm)
plt.loglog(freq, psd_cl / psd_noise_slow)
# plt.loglog(freq, psd_cl / psd_ncp_slow)
plt.xlabel("Frequency (Hz)")
plt.ylabel(r"|ETF|${}^2$");
# %%
