using multiwfs
using QuadGK


pref = target^2 / quadgk(psd_von_karman, 0, 500)[1]
sqrt(quadgk(f -> psd_von_karman(f, pref), 0, 500)[1])

