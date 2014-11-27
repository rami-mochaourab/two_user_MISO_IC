%% calculates the Pareto optimal parameters in closed form: Given parameter lambda2, calculate lambda1
function [la1, la2] = calc_contract_sqrt(H11,H12,H21,H22,n, lambdaMRT1,lambdaMRT2)

g11 = norm(H11)^2;
g12 = norm(H12)^2;
g21 = norm(H21)^2;
g22 = norm(H22)^2;
tg12 = abs(H11'*H12)^2;
tg21 = abs(H22'*H21)^2;
tg11 = g11 - tg12/g12;
tg22 = g22 - tg21/g21;

la2 = 10^-15:lambdaMRT2/2000:lambdaMRT2;

z1 = zeros(1,length(la2));

for idx = 1:1:length(la2)
    
    B = (g21*g12*( sqrt(la2(idx)*(g22-tg22)) + sqrt((1-la2(idx))*tg22) ))/...
        ( ( sqrt((g22-tg22))/sqrt(la2(idx)) - sqrt(tg22)/sqrt(1-la2(idx)) ) * (n+la2(idx)*g21));
    a = -(B-g12)^2 * g11;
    b = (B-g12)*((tg11+g11)*B + (2*n-g12)*(g11-tg11)+2*tg11*n);
    c = (-tg11*B^2 - 2*n*g11*B + n*(2*g12 - n)*(g11-tg11) - tg11*n^2);
    d = (g11-tg11)*n^2;
    
    sols = roots([a, b, c, d]);
    
    for idx2 = 1:1:length(sols)
        if isreal(sols(idx2)) && sols(idx2) >=0
            check = isPareto(H11,H12,H21,H22,n, sols(idx2),la2(idx));
            checks(idx2) = check;
        else
            checks(idx2) = 0;
        end
    end
    
    t =  sols(checks == 1);
    %if length(t) > 1
    %    t
    %end
    z1(idx) = t(1);
    
end

la1 = z1;