%% plot optimal operating points in the two-user MISO Interference channel 
% calculate the Pareto boundary of the rate region and the Walrasian
% equilibrium using the results in the paper:
%
% R. Mochaourab and E. A. Jorswieck "Exchange Economy in Two-User Multiple-Input
% Single-Output Interference Channels" IEEE J. Sel. Topics Signal Process., 
% Special Issue on Game Theory in Signal Processing, vol. 6, no. 2, pp. 151-164, 
% Apr. 2012. 
%
% Email: rami.mochaourab@ee.kth.se
% 
% Date: November 2014
% --------------------------------------------------------------

clear;
clc;
% close all;

%% Initialize parameters

numAntTx = 2;       % number of transmit antennas

noisePower_db = 0;  % variance of the noise in dB
noisePower = 10.^(noisePower_db./10);

rho = 1/noisePower;
rho_db = 10*log10(rho); % signal to noise ration in dB

%% Generate Channels
H = 1/sqrt(2) * (randn(numAntTx, 4) + 1i * randn(numAntTx,4));

% Hjk : channel from Tx j to Rx k
H11 = H(:,1); 
H12 = H(:,3); 
H21 = H(:,4); 
H22 = H(:,2);

%% Paremetrization of optimal Beamforming vectors

samples = 100; % number of parameter samples

lambdaMRT1 = abs(H11'*H12)^2/(norm(H11)*norm(H12))^2;
lambdaMRT2 = abs(H22'*H21)^2/(norm(H22)*norm(H21))^2;

lambda_vect1 = 0:lambdaMRT1/samples:lambdaMRT1;
lambda_vect2 = 0:lambdaMRT2/samples:lambdaMRT2;

[w1_opt] = calc_beams(H11,H12,lambda_vect1);
[w2_opt] = calc_beams(H22,H21,lambda_vect2);

% Calculate Rate tuples which include the Pareto optimal points
dim = length(w2_opt(1,:));

x11 = (H11' * w1_opt) .* (H11' * w1_opt).'';
x12 = (H12' * w1_opt) .* (H12' * w1_opt).'';
x21 = (H21' * w2_opt) .* (H21' * w2_opt).'';
x22 = (H22' * w2_opt) .* (H22' * w2_opt).'';

x11_M = repmat(x11.',1,dim);
x12_M = repmat(x12.',1,dim);

x21_M = repmat(x21,dim,1);
x22_M = repmat(x22,dim,1);

R_pareto1 = log2(1 + x11_M ./ (noisePower + x21_M));
R_pareto2 = log2(1 + x22_M ./ (noisePower + x12_M));


%% calculate the Contract Curve : these are the parameters that achieve the Pareto optimal points

[lambda1c, lambda2c] = calc_contract_sqrt(H11,H12,H21,H22, noisePower,lambdaMRT1,lambdaMRT2);
[lambda2d, lambda1d] = calc_contract_sqrt(H22,H21,H12,H11, noisePower,lambdaMRT2,lambdaMRT1);

% Calculate Pareto boundary

R1_Pareto = zeros(1,length(lambda1c));
R2_Pareto = zeros(1,length(lambda1c));
R1_Pareto2 = zeros(1,length(lambda1c));
R2_Pareto2 = zeros(1,length(lambda1c));

for idx = 1:1:length(lambda1c)
    
    [R1_Pareto(idx), R2_Pareto(idx)] = calc_Rates_lambdas(lambda1c(idx),lambda2c(idx), H11,H12,H21,H22,noisePower);
    [R1_Pareto2(idx), R2_Pareto2(idx)] = calc_Rates_lambdas(lambda1d(idx),lambda2d(idx), H11,H12,H21,H22,noisePower);
    
end

%% Calculate Walrasian equilibrium
[Rate_Walras1, Rate_Walras2] = calc_Walras(H11,H12,H21,H22,noisePower);

%% Plots

RateRegion = figure;
axes('Parent',RateRegion,'FontSize',9,'FontName','times');

title('');
xlabel('R_1(w_1,w_2) [bpcu]','FontSize',9,'FontName','times');
ylabel('R_2(w_1,w_2) [bpcu]','FontSize',9,'FontName','times');

box on;
grid on;
hold on;

figure(RateRegion);
% parametrization
plot(R_pareto1,R_pareto2,'.', 'Color',[0.7 0.7 0.7]);

% strong Pareto boundary
plot(R1_Pareto(1:end),R2_Pareto(1:end), '-','linewidth',1, 'Color', [0 0 0]);
plot(R1_Pareto2(1:end),R2_Pareto2(1:end), '-','linewidth',1, 'Color', [0 0 0]);

% Weak Pareto boundary
line([R1_Pareto(1) R1_Pareto(1)],[0 R2_Pareto(1)],'linewidth',1, 'Color', [0 0 0])
line([0 R1_Pareto2(1)],[R2_Pareto2(1) R2_Pareto2(1)],'linewidth',1, 'Color', [0 0 0])

% Walrasian Equilibrium
plot(Rate_Walras1,Rate_Walras2, '+','linewidth',2,'MarkerSize', 8, 'Color', [0.6 0 0]);
text(Rate_Walras1,Rate_Walras2,'   Walrasian equilibrium', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')

% Nash equilibrium
plot(R_pareto1(end,end), R_pareto2(end,end), 'x','linewidth',2 ,'MarkerSize', 8, 'Color', [0 0 0]);
text(R_pareto1(end,end), R_pareto2(end,end),'   Nash equilibrium', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')

% joint zero forcing
plot(R_pareto1(1,1), R_pareto2(1,1), 'x','linewidth',2,'MarkerSize', 8, 'Color', [0 0 0]);
text(R_pareto1(1,1), R_pareto2(1,1),'   joint zero forcing', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
