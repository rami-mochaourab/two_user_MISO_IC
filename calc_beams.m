%% calculate the beamforming vector w for given parameter lambda
% H_direct : direct channel
% H_interference : interference channel
%
function [w] = calc_beams(H_direct, H_interference, lambda_v)

w = zeros(length(H_direct),length(lambda_v));

Proj1 = ((H_interference * H_interference') / (H_interference'*H_interference));
Proj2 = ( eye(length(H_interference)) - (H_interference * H_interference') / (H_interference'*H_interference));

for idx_lambda = 1:1:length(lambda_v)
    lambda = lambda_v(idx_lambda);
    w(:,idx_lambda) = sqrt(lambda) * Proj1 * H_direct / norm(Proj1 * H_direct)  + sqrt(1-lambda) * (Proj2 * H_direct) / norm(Proj2 * H_direct);
end