import numpy as np
from scipy.signal import windows, welch

def rms(data, axis=0, places=8):
	"""
	Computes the root-mean-square of `data` to `places` places.
	"""
	return round(
		np.sqrt( # root
			np.mean( # mean
				(data - np.mean(data, axis=axis)) ** 2 # square
			)
		),
		places
	)

def genpsd(tseries, dt, nseg=4, remove_dc=True):
	nperseg = 2**int(np.log2(tseries.shape[0]/nseg))
	# firstly ensures that nperseg is a power of 2
	# secondly ensures that there are at least nseg segments 
	# per total time series length for noise averaging
	window = windows.hann(nperseg)
	freq, psd = welch(
		tseries,
		fs=1./dt,
		window=window,
		noverlap=nperseg*0.25,
		nperseg=nperseg,
		detrend=False,
		scaling='density'
	)
	if remove_dc:
		freq, psd = freq[1:], psd[1:] #remove DC component (freq=0 Hz)
	return freq, psd