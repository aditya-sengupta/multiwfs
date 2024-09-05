using multiwfs
using Symbolics
using Polynomials
using Plots

f_loop = 1000.0

freq_low = 0.0079
damp_low = 1.0
freq_high = 2.7028
damp_high = 0.2533
log_lf_cost = 1.7217
log_lf_process_noise = -8.0
log_hf_cost = -0.6807
log_hf_process_noise = -8.0

Av1 = A_vib(freq_high/f_loop, damp_high)
Av2 = A_vib(freq_low/f_loop, damp_low)
A_ar1 = [0.995 0; 1 0]
L = A_DM(2)
Ã = block_diag(L, A_ar1, Av1, Av2)
C̃ = [0 -1 0 1 0 1 0 1]
D̃ = [1 0 0 0 0 0 0 0]' 
B = [0; 0; 1; 0; exp10(log_hf_process_noise); 0; exp10(log_lf_process_noise); 0]
Pw = hcat(1...)
W = B * Pw * B'
V = hcat(1...)
K̃ = kalman_gain(Ã, C̃, W, V)
Vv = [0 -1 0 1 0 exp10(log_hf_cost) 0 exp10(log_lf_cost)]
Q = Vv' * Vv
R = zeros(1,1)
L = lqr_gain(Ã, D̃, Q, R)

#analytic_tf = simplify(lqg_controller_tf(Ã, D̃, C̃, K̃, L, 1/x))

fr = exp10.(-4:0.01:log10(f_loop/2))
sr = 2π .* im * (fr ./ f_loop)
zr = exp.(sr)
#[substitute(analytic_tf, Dict(x => 1/z)) for z in zr]
#[lqg_controller_tf(Ã, D̃, C̃, K̃, L, z) for z in zr]

"""function extract_coeffs(symbolic_args)
	i = 0
	coeffs = []
	while true
		new_coeff = Symbolics.coeff(symbolic_args, x^i)
		if (i > 0) & (new_coeff ≈ 0)
			return coeffs
		end
		push!(coeffs, new_coeff)
		i += 1
	end
end"""

#numerator_coeffs = Float64.(extract_coeffs(simplify(analytic_tf.val.arguments[1], expand=true)))
#denominator_coeffs = Float64.(extract_coeffs(simplify(analytic_tf.val.arguments[2], expand=true)))
#normalizer = numerator_coeffs[end]
#numerator_coeffs ./= normalizer
#denominator_coeffs ./= normalizer
#lqg_zeros = roots(Polynomial(numerator_coeffs))
#lqg_poles = roots(Polynomial(denominator_coeffs))

search_gain!(sys_high)
nyquist_plot(sys_high)

plot(fr, abs2.(1 ./ (1 .+ lqg_ol_tf)), xscale=:log10, yscale=:log10)