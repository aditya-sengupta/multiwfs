using multiwfs
using DSP

white_noise = rand(Normal(), 1000000);
begin
    f_loop = 1000
    8 * sum(power(psd(white_noise, f_loop))) * f_loop / length(white_noise), var(white_noise)
end

