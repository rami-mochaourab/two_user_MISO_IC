%%
function Bin_rel = isPareto(~,H12,H21,H22,n, lambda1,lambda2)

g12 = norm(H12)^2;
g21 = norm(H21)^2;
g22 = norm(H22)^2;

tg21 = abs(H22'*H21)^2;
tg22 = g22 - tg21/g21;

la2 = lambda2;
la1 = lambda1;

C = (( sqrt(la2*(g22-tg22)) + sqrt((1-la2)*tg22) ))/...
        ( ( sqrt((g22-tg22))/sqrt(la2) - sqrt(tg22)/sqrt(1-la2) ) * (n/g21+la2));

cond1 = (n/g12 + la1 -C*la1 <= 0 && n/g12 + la1 -C*la1 +C <= 0);
cond2 = (n/g12 + la1 -C*la1 >= 0 && n/g12 + la1 -C*la1 +C >= 0);

if (cond1 || cond2)
    Bin_rel = true;
else
    Bin_rel = false;
end