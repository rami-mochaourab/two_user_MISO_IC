%% calculater user demands given the prices
function [x11,x12,x21,x22] = calc_Demands(H11,H12,H21,H22,noisePower, p)

g11 = norm(H11)^2;
g12 = norm(H12)^2;
g21 = norm(H21)^2;
g22 = norm(H22)^2;
tg12 = abs(H11'*H12)^2;
tg21 = abs(H22'*H21)^2;
tg11 = g11 - tg12/g12;
tg22 = g22 - tg21/g21;

lambda1NE = (g11 - tg11)/g11;
lambda2NE = (g22 - tg22)/g22;

B = noisePower/g21 + lambda2NE - lambda1NE * p;
x11 = (1)/(1 + tg11/(g11-tg11) * (p/B + 1)^2);
% x12 = lambda2NE - lambda1NE * p + x11 * p; 
x12 = lambda1NE * p - x11 * p; 

C = noisePower/g12 + lambda1NE - lambda2NE /p;
x22 = (1)/(1 + tg22/(g22-tg22) * (1/(p*C) + 1)^2);
% x21 = lambda1NE - lambda2NE / p + x22 /p; 
x21 = lambda2NE / p - x22 /p; 