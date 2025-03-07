function make_double_integrator_system(α, β; f_loop=100.0)
    Vc = 1
    A_ar1 = [α β; 1 0]
    L = A_DM(2)
    Ã = block_diag(L, A_ar1)
    C̃ = [0 -1 0 1]
    D̃ = [1 0 0 0]' 
    B = [0; 0; 1; 0]
    Pw = hcat(1...)
    W = B * Pw * B'
    V = hcat(Vc...)
    K̃ = kalman_gain(Ã, C̃, W, V)
    Vv = [0 -1 0 1]
    Q = Vv' * Vv
    R = zeros(1,1)
    L = lqr_gain(Ã, D̃, Q, R)
    lqg = LQG(Ã, D̃, C̃, K̃, L)
    return AOSystem(f_loop, 0.1, 1.77, 0.999, 1, lqg)
end

export make_double_integrator_system