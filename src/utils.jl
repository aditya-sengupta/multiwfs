using SciPy: signal
using FFTW

function block_diag(matrices...)
    matrices = collect(matrices)
    total_dim = 0
    for m in matrices
        @assert ndims(m) == 2
        @assert size(m, 1) == size(m, 2)
        total_dim += size(m, 1)
    end
    types = [typeof(m).parameters[1] for m in matrices]
    if ComplexF64 in types
        t = ComplexF64
    else
        t = Float64
    end
    A = zeros(t, total_dim, total_dim)
    i = 1
    for m in matrices
        k = size(m, 1)
        A[i:(i+k-1), i:(i+k-1)] .= m
        i += k
    end
    A
end

function genpsd(tseries, f_loop, nseg=4)
	nperseg = 2^Int(round(log2(length(tseries)/nseg))) #firstly ensures that nperseg is a power of 2, secondly ensures that there are at least nseg segments per total time series length for noise averaging
	window = signal.windows.hann(nperseg)
	freq, psd = signal.welch(tseries, fs=f_loop,window=window, noverlap=nperseg*0.25,nperseg=nperseg, detrend=false,scaling="density")
	freq, psd = freq[2:end], psd[2:end] #remove DC component (freq=0 Hz)
	return freq, psd
end

function gen_avg_per_unb(sig1d0, per_len)
    """
    Generate the Averaged, Unbiased Periodogram for a 1D signal using 
    segment length per_len. Options that affect quality include mean removal, 
    windowing, and half-overlapping to improve SNR.

    Args:
        sig1d0 (Vector{Float64}): the 1D signal that we want to estimate the PSD of
        per_len (Int): the segment length to calculate the FFT on

    Returns:
        Vector{Float64}: the 1D periodogram of length per_len
    """
    total_len = length(sig1d0)
    if per_len > total_len
        throw(ArgumentError("Error! Input periodogram length $per_len is longer than data length $total_len!"))
    end

    sig1d = sig1d0 .- mean(sig1d0)
    num_intervals = floor(Int, total_len / (per_len / 2)) - 1
    start_indices = collect(0:(num_intervals - 1)) .* div(per_len, 2)
    ind = 0:(per_len - 1)
    window = 0.5 .- 0.5 .* cos.(2π .* ind ./ (per_len - 1))
    psd = zeros(Float64, per_len)

    for a in 1:num_intervals
        segment = sig1d[start_indices[a] + 1 : start_indices[a] + per_len]
        psd .+= abs.(FFTW.fft(window .* segment)).^2 ./ per_len^2
    end

    psd ./= num_intervals
    psd ./= sum(window.^2) / per_len
    return psd
end

no_filter = ZPKFilter(0, 0, 1)

subtract_mean(x) = x .- mean(x)
rms(x) = sqrt(mean(x .^ 2))

export block_diag, genpsd, gen_avg_per_unb, no_filter, subtract_mean, rms