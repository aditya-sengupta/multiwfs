using LinearAlgebra: I
using ControlSystems: are, Discrete
using SciPy: linalg

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
    return (L * inv(exp(s)*I - A + B*L + K*C) * K)[1,1]
end

function kalman_gain(A, C, W, V)
    P = are(Discrete(1), A', C', W, V)
    K = P * C' * inv(C * P * C' + V)
    return K
end

function lqr_gain(A, B, Q, R)
    P = nothing
    try
        P = are(Discrete(1), A, B, Q, R)
    catch 
        P = linalg.solve_discrete_are(A, B, Q, R)
    end
    L = -inv(B' * P * B) * B' * P * A
    return L
end

export A_vib, dynamic_system_tf, lqg_tf, kalman_gain, lqr_gain