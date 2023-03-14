a0 = - 0.5;
b0 = 1;
ND = 100;
y = zeros(ND,1);
y(1,1) = 1;
u  =  randn(ND,1); % buzení systému
e  =  randn(ND,1); % bílý šum
e2 =  zeros(ND,1); % barevný šum
e2(1,1) = randn(1,1);
c0 = 0.8;

for i = 2:ND
    e(i,1) = c0*e2(i-1, 1);
    y(i,1) = -a0*y(i-1,1) + b0*u(i,1) + e(i, 1) + e2(i,1);
end

theta = [y(1:ND-1,1), u(1:ND-1,1)];
E = e(1:ND-1, 1);

odhad = inv(theta'* theta)*theta'.* y(ND-1,1)';

theta_odhad = theta' + inv(theta'* theta)*theta'* e(ND-1,1)';
chyba = theta'  - theta_odhad;


mse(chyba(1,:))

mse(chyba(2,:))


