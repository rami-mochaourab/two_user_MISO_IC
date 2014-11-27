function SINR = calc_SINR(H11,H21,w1_d,w2_d, noisePower)

SINR = abs(H11'*w1_d)^2/(noisePower + abs(H21'*w2_d)^2);