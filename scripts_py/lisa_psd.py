import numpy as np

def genAvgPerUnb(sig1d0, per_len, meanrem=True, hanning=True, halfover=True):
    """
    Based on IDL's gen_avg_per_unb.pro.
    Generate the Averaged, Unbiased Periodogram for a 1D signal using 
    segment length per-len. OPtions that affect quality include mean removal, 
    windowing and half-overlapping to improve SNR.

    Args:
        sig1d (float): the 1D signal that we want to estimate the PSD of
        per_len (int): the segment length to calculate the FFT on
        meanrem (bool): If True (default), remove the average value over the entire sig1d first. Improve quality near DC.
        hanning (bool): If True (default), use ahanning window to improve dynamic range (at the cost of resolution)
        halfover (bool): If True (default), use half-overlapped segments to improve SNR on estimate

    Returns:
        the 1D periodogram of length per_len.

    Example:
        ::

            # need an example here!

    """
    total_len = len(sig1d0)
    if per_len > total_len:
        raise Exception('Error! input periodogram length %i is longer than data length %i!' %(per_len, total_len))

    if meanrem:
        sig1d = sig1d0 - np.sum(sig1d0)/total_len
    else:
        sig1d = sig1d0

    if halfover:
        num_intervals = np.floor(total_len/(per_len/2)) - 1
        start_indices = np.arange(0, int(num_intervals))*int(per_len/2)
    else:
        num_intervals = np.floor(total_len/(per_len)) 
        start_indices = np.arange(0, int(num_intervals))*int(per_len)

    ind = np.arange(0, int(per_len))
    if hanning:
        window = 0.5 - 0.5*np.cos(2*np.pi*ind/(per_len-1))
    else:
        window = np.ones((per_len))

    psd = np.zeros((per_len))
    for a in range(int(num_intervals)):
        psd = psd + (np.abs(np.fft.fft(window*sig1d[start_indices[a]:(start_indices[a] + per_len - 1 + 1)]))**2)/per_len**2

    psd = psd/num_intervals
    psd = psd/(np.sum(window**2)/per_len)
    return psd


def genAvgPerUnb(sig1d0, per_len):
    """
    Based on IDL's gen_avg_per_unb.pro.
    Generate the Averaged, Unbiased Periodogram for a 1D signal using 
    segment length per-len. OPtions that affect quality include mean removal, 
    windowing and half-overlapping to improve SNR.

    Args:
        sig1d (float): the 1D signal that we want to estimate the PSD of
        per_len (int): the segment length to calculate the FFT on
        meanrem (bool): If True (default), remove the average value over the entire sig1d first. Improve quality near DC.
        hanning (bool): If True (default), use ahanning window to improve dynamic range (at the cost of resolution)
        halfover (bool): If True (default), use half-overlapped segments to improve SNR on estimate

    Returns:
        the 1D periodogram of length per_len.

    Example:
        ::

            # need an example here!

    """
    total_len = len(sig1d0)
    if per_len > total_len:
        raise Exception('Error! input periodogram length %i is longer than data length %i!' %(per_len, total_len))

    sig1d = sig1d0 - np.sum(sig1d0)/total_len
    num_intervals = np.floor(total_len/(per_len/2)) - 1
    start_indices = np.arange(0, int(num_intervals))*int(per_len/2)
    ind = np.arange(0, int(per_len))
    window = 0.5 - 0.5*np.cos(2*np.pi*ind/(per_len-1))
    psd = np.zeros((per_len))
    for a in range(int(num_intervals)):
        psd = psd + (np.abs(np.fft.fft(window*sig1d[start_indices[a]:(start_indices[a] + per_len - 1 + 1)]))**2)/per_len**2

    psd = psd/num_intervals
    psd = psd/(np.sum(window**2)/per_len)
    return psd
