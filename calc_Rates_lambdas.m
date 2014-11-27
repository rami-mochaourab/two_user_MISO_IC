%% calculate rates for given parameters lambda1, lambda2
function [R1,R2] = calc_Rates_lambdas(lambda1, lambda2, H11, H12, H21, H22, noisePower)

w1 = calc_beams(H11, H12, lambda1);
w2 = calc_beams(H22, H21, lambda2);

R1 = log2(1+abs(H11'*w1)^2 / (noisePower + abs(H21'*w2)^2));
R2 = log2(1+abs(H22'*w2)^2 / (noisePower + abs(H12'*w1)^2));