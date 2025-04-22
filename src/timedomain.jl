"""
Generate a time series with a given power spectral density (PSD) model.

# Arguments
- `psd_model::Function`: A function that takes a frequency and returns the PSD value at that frequency.
- `fmax::Float64`: The maximum frequency of interest.
- `npoints::Int`: The number of points in the discrete PSD.

# Returns
- `time_series::Vector{Float64}`: The generated time series.

copilot-generated off this:
https://dsp.stackexchange.com/questions/76660/generate-a-time-series-from-power-spectral-density 

This doesn't seem to produce the correct power for white noise...
"""
function generate_time_series(psd_model::Function, fmax::Float64, npoints::Int)
    fs = 2 * fmax
    df = fs / npoints
    freqs = df * (0:(npoints ÷ 2))
    psd_values = psd_model.(freqs)
    amplitudes = sqrt.(npoints * fs * psd_values)
    phases = 2π * rand(length(amplitudes))
    complex_amplitudes = amplitudes .* exp.(-im .* phases)
    full_spectrum = vcat(complex_amplitudes, conj(reverse(complex_amplitudes[2:end-1])))
    time_series_complex = ifft(full_spectrum)
    time_series = real(time_series_complex)
    return time_series
end

function generate_time_series(vk::VonKarman, fmax::Float64, npoints::Int)
    return generate_time_series(f -> psd_von_karman(f, vk), fmax, npoints)
end

function generate_openloop_timeseries(sim, N)
    atm_timeseries = generate_time_series(sim.vk_atm, sim.f_loop/2, N)
	slow_ncp_timeseries = subtract_mean(generate_time_series(sim.vk_ncp, sim.f_loop/2, N))
	fast_ncp_timeseries = subtract_mean(generate_time_series(sim.vk_ncp, sim.f_loop/2, N))
	fast_noise_timeseries = subtract_mean(generate_time_series(f -> psd_von_karman(sim.f_noise_crossover, sim.vk_atm), sim.f_loop/2, N))
	slow_noise_timeseries = subtract_mean(generate_time_series(f -> psd_von_karman(sim.f_noise_crossover, sim.vk_atm), sim.f_loop/2, N))
    return (atm=atm_timeseries, slowncp=slow_ncp_timeseries, fastncp=fast_ncp_timeseries, slownoise=slow_noise_timeseries, fastnoise=fast_noise_timeseries)
end

function timedomain_closedloop(sim, all_timeseries)
	N = length(all_timeseries.atm)
    fastcon, slowcon = sim.fast_controller, sim.slow_controller
    phase_seen_by_slow_sensor = zeros(N)
	global fast_last_c = 0.0
	global slow_last_c = 0.0
	global slow_this_s = 0.0
	global slow_this_c = 0.0
	global current_dm = 0.0
	fast_buffer = [0.0, 0.0]
	slow_buffer = zeros(1+sim.R)
	for i in 1:N
		slow_update_this_iteration = i % sim.R == 0
		this_e = all_timeseries.atm[i] - current_dm
		fast_buffer[end] = this_e + all_timeseries.fastncp[i]
		fast_this_s = fast_buffer[begin]
		fast_this_s += all_timeseries.fastnoise[i]
		fast_this_s = output!(fastcon.cfilter, fast_this_s)[1]
		fast_buffer[1:end-1] .= fast_buffer[2:end]
		slow_buffer[end] = this_e + all_timeseries.slowncp[i]
        phase_seen_by_slow_sensor[i] = slow_buffer[begin]
		if slow_update_this_iteration
			global slow_this_s = sum(slow_buffer[1:end-1])/sim.R + all_timeseries.slownoise[i]
		end
		slow_buffer[1:end-1] = slow_buffer[2:end]
		fast_this_c = fastcon.gain * fast_this_s + fastcon.leak * fast_last_c
		if slow_update_this_iteration
			global slow_this_c = slowcon.gain * slow_this_s + slowcon.leak * slow_last_c
		end
		global fast_last_c = fast_this_c
		global slow_last_c = slow_this_c
		global current_dm = fast_this_c + slow_this_c
	end
	return phase_seen_by_slow_sensor
end

export generate_time_series, generate_openloop_timeseries, timedomain_closedloop