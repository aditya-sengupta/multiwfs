using LinearAlgebra: I

function A_vib(f_over_f_loop, k)
    ω = 2π * f_over_f_loop
    a₁ = 2 * exp(-k * ω) * cos(ω * sqrt(Complex(1 - k^2)))
    a₂ = -exp(-2 * k * ω)
    [a₁ a₂; 1 0]
end

function dynamic_system_tf(s, A, B, C)
    return (C * inv(exp(s)*I - A) * B)[1,1]
end

function lqg_tf(s, A, B, C, K, L)
    return (L * inv(s*I - A + B*L + K*C) * K)[1,1]
end

export A_vib, dynamic_system_tf, lqg_tf