import numpy as np
from scipy.signal import welch, windows
from matplotlib import pyplot as plt
from matplotlib import cm

def f2s(f):
	return 1j * 2.0 * np.pi * f

def genpsd(tseries, dt, nseg=4):
	nperseg = 2**int(np.log2(tseries.shape[0]/nseg)) #firstly ensures that nperseg is a power of 2, secondly ensures that there are at least nseg segments per total time series length for noise averaging
	window = windows.hann(nperseg)
	freq, psd = welch(tseries, fs=1./dt,window=window, noverlap=nperseg*0.25,nperseg=nperseg, detrend=False,scaling='density')
	freq, psd = freq[1:], psd[1:] #remove DC component (freq=0 Hz)
	return freq, psd

def get_freq(f_loop, N):
	min_freq = f_loop / (2 * N)
	return np.arange(min_freq, f_loop / 2 + min_freq, min_freq)

class AOSystem:
	def __init__(self, f_loop, frame_delay, gain, leak, fpf=10, filter_type="none", filter_cutoff=None):
		self.f_loop = f_loop
		self.Ts = 1 / f_loop
		self.frame_delay = round((frame_delay-np.floor(frame_delay))*fpf)/fpf+np.floor(frame_delay)
		self.tau = self.frame_delay * self.Ts
		self.gain = gain
		self.leak = leak
		self.filter_type = filter_type
		if self.filter_type != "none":
			self.alpha = np.exp(-2*np.pi*filter_cutoff/f_loop)
		
	def Hwfs(self, s):
		return (1.0 - np.exp(-self.Ts * s)) / (self.Ts * s)

	def Hzoh(self, s):
		return self.Hwfs(s)

	def Hlag(self, s):
		return np.exp(-self.tau * s)

	def Hint(self, s):
		return self.gain / (1.0 - np.exp(-self.Ts * s)) #integrator

	def Hlint(self, s):
		return self.gain / (1.0 - self.leak * np.exp(-self.Ts * s)) #leaky integrator

	def Hcont(self, s):
		return self.Hlint(s)

	def Hfilter(self, s):
		if self.filter_type == "high":
			return self.alpha * (1 - np.exp(-s / self.f_loop)) / (1 - self.alpha * np.exp(-s / self.f_loop))
		elif self.filter_type == "low":
			return (1 - self.alpha) / (1 - self.alpha * np.exp(-s / self.f_loop))
		else:
			return 1

	def Holsplane(self, s):
		return self.Hwfs(s) * self.Hlag(s) * self.Hcont(s) * self.Hfilter(s) * self.Hzoh(s)

	def Hol(self, f):
		assert np.dtype(f[0]) != "complex128"
		return self.Holsplane(f2s(f)) #open loop transfer function

	def Hrej(self, f):
		return 1.0 / (1.0 + self.Hol(f)) #rejection transfer function

	def Hcl(self, f):
		Hol_val = self.Hol(f)
		return Hol_val / (1.0 + Hol_val)

	def Hn(self, f):
		return self.Hcl(f) / self.Hwfs(f2s(f))

	def plot_frequency_response(self, f):
		F_ol, F_cl = np.abs(self.Hol(f)) ** 2, np.abs(self.Hcl(f)) ** 2
		fig, ax = plt.subplots()
		ax.plot(f, F_ol, label="Open-loop PSD")
		ax.plot(f, F_cl, label="Closed-loop PSD")
		ax.set_xlabel("Temporal frequency (Hz)")
		ax.set_ylabel("PSD (arbitrary units)")
		ax.set_xscale("log")
		ax.set_yscale("log")
		ax.legend()
		axr = plt.gca().twinx()
		axr.set_xscale("log")
		axr.set_yscale("log")
		axr.set_ylabel(r"$\text{rejection transfer function (RTF)}|^2$")
		axr.plot(f, F_cl / F_ol, color="k", label=r"$|\text{RTF}|^2$")
		fig.suptitle(f"Ts = {self.Ts*1e3}ms, tau = {self.tau*1e3}ms, gain = {self.gain}, leak = {self.leak}", size=15)

		axr.legend()
		
		plt.show()
  
	def nyquist(self, f):
		linfreq = np.linspace(np.min(f), np.max(f), 1000)
		self.nyquist_p = self.Hol(linfreq)
		self.nyquist_n = self.Hol(-1*linfreq)
		return linfreq, self.nyquist_p, self.nyquist_n

	def phase_margin_and_location(self):
		phase_margin_point = np.abs(np.real(self.nyquist_p)**2 + np.imag(self.nyquist_p)**2 - 1)
		phase_margin_ind = np.argmin(phase_margin_point)
		return float(round(np.angle(self.nyquist_p[phase_margin_ind]) % np.pi * (180 / np.pi), 2)), self.nyquist_p[phase_margin_ind]

	def gain_margin_and_location(self):
		gain_margin_point = np.abs(np.imag(self.nyquist_p))
		gain_margin_ind = np.argmin(gain_margin_point)
		return float(round(-1 / np.real(self.nyquist_p[gain_margin_ind]), 2)), self.nyquist_p[gain_margin_ind]
  
	def plot_nyquist(self, f, gain_margin_limit=2.5):
		linfreq, Hol2plotp, Hol2plotn = self.nyquist(f)
		colors = cm.viridis(np.linspace(0,0.75,4))
		fig, ax = plt.subplots(figsize=(5,5))
		ax.grid('on')
		ax.axvline(0, color='k', ls='--')
		ax.axhline(0, color='k', ls='--')
		ax.plot(np.real(Hol2plotp), np.imag(Hol2plotp), color=colors[1])
		ax.plot(np.real(Hol2plotn), np.imag(Hol2plotn), color=colors[1])
		phasegrid = np.linspace(-np.pi,np.pi,500)
		xunit, yunit = np.cos(phasegrid), np.sin(phasegrid)
		ax.plot(xunit, yunit, color=colors[2], ls=':')
		ax.axvline(-1/gain_margin_limit,color=colors[0],ls='--')
		ax.plot(-1,0,'-ro')
		ax.set_xlim(-1.1, 1.1)
		ax.set_ylim(-1.1, 1.1)
		ax.plot(np.linspace(-2,0,10), np.linspace(-2,0,10), color=colors[3], ls='--')
		ax.plot(np.linspace(-2,0,10), np.linspace(2,0,10), color=colors[3], ls='--')

		pm, pmloc = self.phase_margin_and_location()
		gm, gmloc = self.gain_margin_and_location()
		ax.plot(np.real(pmloc), np.imag(pmloc), '-o', color=colors[3])
		ax.set_title(f"{gm=}, {pm=}")
		plt.show()
  