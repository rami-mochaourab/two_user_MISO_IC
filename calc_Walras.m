%% calculater the rates in Walrasian equilibrium
function [R1,R2] = calc_Walras(H11,H12,H21,H22,noise)

g11 = norm(H11)^2;
g12 = norm(H12)^2;
g21 = norm(H21)^2;
g22 = norm(H22)^2;
tg12 = abs(H11'*H12)^2;
tg21 = abs(H22'*H21)^2;
tg11 = g11 - tg12/g12;
tg22 = g22 - tg21/g21;

lambdaMRT1 = abs(H11'*H12)^2/(norm(H11)*norm(H12))^2;
lambdaMRT2 = abs(H22'*H21)^2/(norm(H22)*norm(H21))^2;

%% Calculate Walras Roots

n = noise;
g1 = g11-tg11;
g2 = g22-tg22;

a1 = ((g1-tg11)/(g1+tg11)) * (n/g12+lambdaMRT1)^2 * (1-lambdaMRT1)*lambdaMRT1;


a2 = - 2*(1-lambdaMRT1)*lambdaMRT1*(lambdaMRT2+n/g21)*(n/g12+lambdaMRT1)^2 ...
     - 2*(1-lambdaMRT1)*lambdaMRT1*((g1-tg11)/(g1+tg11))*((g2-tg22)/(g2+tg22))*(n/g12+lambdaMRT1);
   
   
a3 = + 2*(tg11/g11-2*tg11/g11*lambdaMRT1+lambdaMRT1^2)*(n/g12+lambdaMRT1)*(1-lambdaMRT2)*lambdaMRT2...
     + 4*((g2-tg22)/(g2+tg22))*(lambdaMRT2+n/g21)*(n/g12+lambdaMRT1)*(1-lambdaMRT1)*lambdaMRT1...
     +   ((g1-tg11)/(g1+tg11))*(lambdaMRT2^2-2*tg22/g22*lambdaMRT2+tg22/g22)*(1-lambdaMRT1)*lambdaMRT1;


a4 = -4*((g1-tg11)/(g1+tg11))*(lambdaMRT2+n/g21)*(n/g12+lambdaMRT1)*(1-lambdaMRT2)*lambdaMRT2...
      - (tg11/g11-2*tg11/g11*lambdaMRT1+lambdaMRT1^2)*((g2-tg22)/(g2+tg22))*(1-lambdaMRT2)*lambdaMRT2...
      - 2*(lambdaMRT2+n/g21)*(lambdaMRT2^2-2*tg22/g22*lambdaMRT2+tg22/g22)*(1-lambdaMRT1)*lambdaMRT1;
  
a5 = + 2*(lambdaMRT2+n/g21)^2*(n/g12+lambdaMRT1)*(1-lambdaMRT2)*lambdaMRT2...
     + 2*((g1-tg11)/(g1+tg11))*((g2-tg22)/(g2+tg22))*(lambdaMRT2+n/g21)*(1-lambdaMRT2)*lambdaMRT2;


a6 = -((g2-tg22)/(g2+tg22)) * (lambdaMRT2+n/g21)^2 * (1-lambdaMRT2)*lambdaMRT2;

roots_Prices_Walras = roots([a1,a2,a3,a4,a5,a6]);

% Conditions on p
ubound_prices = (noise/g21 + lambdaMRT2)/(lambdaMRT1);
lbound_prices = (lambdaMRT2)/(noise/g12 + lambdaMRT1);

WP = roots_Prices_Walras > lbound_prices & roots_Prices_Walras < ubound_prices;
roots_Prices_Walras(WP);

%% Calculate Demands
price_ratio = max(roots_Prices_Walras(WP));
[d11,~,~,d22] = calc_Demands(H11,H12,H21,H22,noise, price_ratio);

w1_d1 = calc_beams(H11, H12, d11);
w2_d2 = calc_beams(H22, H21, d22);
SINR1 = calc_SINR(H11,H21,w1_d1,w2_d2, noise);
SINR2 = calc_SINR(H22,H12,w2_d2,w1_d1, noise);

R1 = log2(1+SINR1);
R2 = log2(1+SINR2);

