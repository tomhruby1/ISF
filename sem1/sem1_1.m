% A
clear; clc; close all;

a11 = -0.5; b11 = 0.9; c11 = 0; lambda = 1;
a12 = -0.5; b12 = 0.9; c12 = -0.4;

ND = 501;
t = linspace(0,ND-1,ND);

problem_size = 3;

y1 = zeros(ND,problem_size);
y1(1,1) = 0;
y2 = zeros(ND,problem_size);
y2(1,1) = 0;

sigma = 1;

e  =  randn(ND,1)*lambda; % bílý šum

u = zeros(ND, 3);
for problem = 1:3
    switch problem
        case 1
            u(:,problem)  =  ones(ND,1); % buzení systému
        case 2
            u(:,problem) = zeros(ND,1);
            u(1,problem) = sigma;
        case 3
            u(:,problem) = randn(ND, 1)*sigma;
    end
    
    for i = 2:ND
        y1(i,problem) = -a11*y1(i-1,problem) + b11*u(i-1,problem) + e(i, 1);
        y2(i,problem) = -a12*y2(i-1,problem) + b12*u(i-1,problem) + e(i, 1) + c12*e(i-1,1);
    end
end

figure;
plot(t, y1);
title('System1');
legend('step', 'impulse', 'gaussian');

figure;
plot(t, y2);
title('System2');
legend('step', 'impulse', 'gaussian');

%% B

theta_odhad1 = zeros(2,3);
theta_odhad2 = zeros(2,3);

for problem = 1:3
    phi = [-y1(1:ND-1,problem), u(1:ND-1,problem)];
    %E = e(1:ND, 1);
    theta_odhad1(:,problem) = inv(phi'* phi)*phi' * y1(2:ND,problem);
    
    phi = [-y2(1:ND-1,problem), u(1:ND-1,problem)];
    theta_odhad2(:,problem) = inv(phi'* phi)*phi' * y2(2:ND,problem);
end

theta_odhad1
theta_odhad2

%% C??

% impuls je blbej, protoze je skoro porad 0, a tim padem asi nic moc nepovi
% o b-koeficientu

% tim padem step je asi dobrej, protoze se dobre projevi A

% gauss se zda byt fakt dobrej (Niki fandi spise stepu!), protoze zkusime ruzne pary u -> y


%% D

y1_reg = zeros(ND,1);
e  =  randn(ND,1)*lambda;


kappa = 0;
v  =  randn(ND,1)*kappa;

k = 2;

% vyber odhadnuteho param
a = theta_odhad1(1,1) 
b = theta_odhad2(1,1)

u = zeros(ND,1)
for i = 2:ND
    u(i) = -k*y1_reg(i-1) + v(i-1);
    y1_reg(i) = -a*y1(i-1) + b*u(i) + e(i, 1);
end


plot(t, y1_reg)

%%
%theta_odhad = phi' + inv(phi'* phi)*phi'* e(ND-1,1)';
chyba = [a11, b11]'  - theta_odhad;