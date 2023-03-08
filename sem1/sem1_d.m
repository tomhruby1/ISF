% D

% - asi generujueme ruzne odhady v zavislosti na ruznych kappa: 
% nebo teoreticky???

clear; clc; close all;

a11 = -0.5; b11 = 0.9; c11 = 0; lambda = 1;

ND = 501;
t = linspace(0,ND-1,ND);

y1 = zeros(ND,1);

e = randn(ND,1)*lambda; % bílý šum

N_kappa = 1000;
max_kappa = 20;


kappas = linspace(0,max_kappa,N_kappa);

theta_odhady = zeros(2, N_kappa);
us = zeros(ND, N_kappa);
for q = 1:N_kappa
    
    kappa = kappas(q);
    
    v  =  randn(ND,1)*kappa;

    k = 0.8; % zafixovano tak, aby pro kappa = 0 stabilni

    u = zeros(ND,1);
    for i = 2:ND
        u(i) = -k*y1(i-1) + v(i-1);
        y1(i) = -a11*y1(i-1) + b11*u(i) + e(i, 1);
    end
    
    us(:,q) = u; 
    
%     figure;
%     plot(t, y1);
%     title('System1 feedback');
    
    phi = [-y1(1:ND-1), u(1:ND-1)];
    theta_odhady(:,q) = inv(phi'* phi)*phi' * y1(2:ND);
   
end

figure; hold on
plot(kappas, theta_odhady(1,:))
plot(kappas, theta_odhady(2,:))
legend('\hat a', '\hat b')

