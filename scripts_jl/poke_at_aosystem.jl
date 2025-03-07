using multiwfs

f_loop = 1000.0
fr = 10 .^ (-2:0.01:log10(f_loop/2))
sr = 2π .* im .* fr
sT = sr / f_loop
no_filter = ZPKFilter(0, 0, 1)
r0 = 0.1031
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
f_noise_crossover = 500.0
noise_normalization = psd_von_karman(f_noise_crossover, vk_atm)
ar1_high = ar1_filter(15.0, f_loop, "high")
sys_fast = AOSystem(f_loop, 1.0, 0.4, 0.999, 1, ar1_high)
sys_slow = AOSystem(f_loop / 10, 0.1, 1.4, 0.999, 1, no_filter)

sys_fast.Ts
sys_fast.frame_delay
sys_fast.τ
sys_slow.fpf

sys_slow.Ts
sys_slow.frame_delay
sys_slow.τ
sys_slow.fpf