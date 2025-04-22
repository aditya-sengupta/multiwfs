# %%
import numpy as np
from matplotlib import pyplot as plt
from lisa_psd import genAvgPerUnb
from ben_lqgic.functions import genpsd
fig, axes = plt.subplots(1, 2, figsize=(9, 6))

for i, f in enumerate(["fastnoise", "fastncp"]):
    atm = np.load(f"../data/olt/olt_{f}.npy")
    kv = np.random.randint(0, 100_000 - 1025, size=(1000,))

    freq_l = np.linspace(0.0, 500.0, 1025)[1:]
    deviations_lisa = []
    deviations_ben = []

    for k in kv:
        ts = atm[k:1025+k]
        rms_truth = np.std(ts)
        sqrt_auc_lisa = float(np.sqrt(np.sum(genAvgPerUnb(ts, 1024))))
        freq_b, psd_b = genpsd(ts, 1e-3)
        sqrt_auc_ben = float(np.sqrt(np.sum(psd_b) * np.diff(freq_b)[0]))
        deviation_lisa = (sqrt_auc_lisa - rms_truth) / rms_truth * 100
        deviation_ben = (sqrt_auc_ben - rms_truth) / rms_truth * 100
        deviations_lisa.append(deviation_lisa)
        deviations_ben.append(deviation_ben)

    mean_deviation_lisa = np.mean(deviations_lisa)
    mean_deviation_ben = np.mean(deviations_ben)

    ax = axes[i]
    ax.hist(deviations_lisa, bins=50, alpha=0.5, label="genAvgPerUnb")
    ax.hist(deviations_ben, bins=50, alpha=0.5, label="genpsd", color="forestgreen")
    ax.axvline(mean_deviation_lisa, color='blue', linestyle='dashed', linewidth=1, label=f"Mean genAvgPerUnb: {mean_deviation_lisa:.2f}%")
    ax.axvline(mean_deviation_ben, color='forestgreen', linestyle='dashed', linewidth=1, label=f"Mean genpsd: {mean_deviation_ben:.2f}%")
    ax.set_title(f"Deviation wrt time-series RMSE, {f}")
    ax.set_ylabel("Count")
    ax.legend()

    ax.set_xlabel("Percent error (%)")
plt.tight_layout()
plt.show()
# %%
k = np.random.randint(0, 100_000 - 1025)
print(k)
ncptimeseries = np.load(f"../data/olt/olt_fastncp.npy")[k:k+2048]
analytic_ncp = 0.25 * (freq_l + 0.001)**(-2)
analytic_noise = (np.ones(len(freq_l)) * 11.02 * (51) ** (-3))

ts = ncptimeseries
a = analytic_ncp
title = "Fast NCP"
freq_b, psd_b = genpsd(ts, 1e-3)
psd_l = genAvgPerUnb(ts, 2048)
sqrt_auc_lisa = float(np.sqrt(np.sum(psd_l)))
sqrt_auc_ben = float(np.sqrt(np.sum(psd_b) * np.diff(freq_b)[0]))
plt.loglog(freq_l, psd_l[:1024], label=f"genAvgPerUnb, sqrt(auc) = {sqrt_auc_lisa:.3f}")
plt.loglog(freq_b, psd_b, label=f"genpsd, sqrt(auc) = {sqrt_auc_ben:.3f}")
plt.loglog(freq_l, a, label="Target PSD", color="black", ls="--")
plt.xlabel("Frequency (Hz)")
plt.ylabel(r"Power (rad${}^2$/Hz)")
plt.title(f"{title}, RMSE = {np.std(ts):.3f}")
plt.legend()

plt.tight_layout()
plt.show()
# %%
a = genAvgPerUnb(ts, 1024)
N = len(a) // 2
np.allclose(a[1:N], np.flip(a[N+1:]))
# %%
