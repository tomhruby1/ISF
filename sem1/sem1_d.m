% A
clear; clc; close all;

a11 = -0.5; b11 = 0.9; c11 = 0; lambda = 1;

ND = 501;
t = linspace(0,ND-1,ND);

y1 = zeros(ND,1);

e = randn(ND,1)*lambda; % bílý šum

N_kappa = 20;
kappas = linspace(0,2,N_kappa);

for i = 1:N_kappa
    
    kappa = kappas(i);
    
    v  =  randn(ND,1)*kappa;

    k = 0.8;


    u = zeros(ND,1);
    for i = 2:ND
        u(i) = -k*y1(i-1) + v(i-1);
        y1(i) = -a11*y1(i-1) + b11*u(i) + e(i, 1);
    end



    figure;
    plot(t, y1);
    title('System1 feedback');

    theta_odhad1 = zeros(2,3);
    theta_odhad2 = zeros(2,3);

    for problem = 1:3
        phi = [-y1(1:ND-1,problem), u(1:ND-1,problem)];
        %E = e(1:ND, 1);
        theta_odhad1(:,problem) = inv(phi'* phi)*phi' * y1(2:ND,problem);
    end

end

theta_odhad1
theta_odhad2

