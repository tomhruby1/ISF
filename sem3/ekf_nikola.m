%% EXTENDED KALMAN FILTER
% inicializace
close all; clear; clc;

dt = 0.02;
N = 455;
ALPHA_Q = 1e-3;
ALPHA_P = 1e-6;

global f B F Q R m_x0 P0

% matice dynamiky
f = [
        1 0 dt 0;
        0 1 0 dt;
        0 0 1 0 ;
        0 0 0 1 ;
    ];
B = [0 0; 0 0; 1 0; 0 1];
% KF mat.
F = [
        1 0 dt 0 0 0;
        0 1 0 dt 0 0;
        0 0 1 0  1 0;
        0 0 0 1  0 1;
        0 0 0 0  1 0;
        0 0 0 0  0 1;
     ]; 
 
Q = 10e-6 * eye(6);
Q(5,5) = ALPHA_Q; Q(6,6) = ALPHA_Q;

R = [
        4e-4 * (pi/180)^2      0;
        0                    1e-4
    ];

m_x0 = [20; 50; 0; -12; 0; 0;];


P0 = eye(6);
P0(5,5) = ALPHA_P; P0(6,6) = ALPHA_P; %pocatecni podminka pro rizeni


%% SIMULATION
% precompute noises
% sim.w = mvgrnd(zeros(size(Q,1),1), Q, N)';
% sim.w_u = mvgrnd([0;0], eye(2)*ALPHA_Q, N)'; 
% sim.v = mvgrnd([0;0], R, N)'; 


sim.w = mvnrnd(zeros(size(Q,1),1), Q, N)';
sim.w_u = mvnrnd([0;0], eye(2)*ALPHA_Q, N)'; 
sim.v = mvnrnd([0;0], R, N)'; 


% RANDOM WALK
% sim.x = m_x0(1:4);
% sim.u = [0; 0];
% sim.z = [0; 0];
% for k = 1:1:N-1
%     [sim.x(:,k+1), sim.u(:,k+1)] = state_eq(sim.x(:,k), sim.u(:,k), sim.w(1:4,k), sim.w_u(:,k)); % dyn. sim.
%      sim.z(:,k+1) = get_observ(sim.x(:,k+1), sim.v(:,k+1)); % measurement sim
%end


% POSKYTNUTA DATA
load('isf_4_data.mat');
sim.x = x;
sim.z = z;
sim.u = u;
N = size(sim.x, 2);

figure(1);
plot(sim.x')
legend('x', 'y', 'x*', 'y*');
title('Stavy');
xlabel('t');
%% estimace


tic
 kfres = kf(N, sim.z);
toc


%% RESULTS
k_start = 1; % plot K-first steps
k_end = 454;
t = k_start:k_end;

figure (2); hold on;
plot(kfres.xp(1,t), kfres.xp(2,t), 'c-');
plot(kfres.xp(1,1), kfres.xp(2,1), 'cx');
plot(kfres.xf(1,t), kfres.xf(2,t), 'b-');
plot(kfres.xf(1,1), kfres.xf(2,1), 'bx');
plot(sim.x(1,t), sim.x(2,t), 'r-');
plot(sim.x(1,1), sim.x(2,1), 'r*');
%legend('pos pred','pos pred', 'pos filt', 'pos filt', 'reálná poloha lodě');
legend('predikovaná pozice','začátek predikce','filtrovaná pozice', 'začátek filtrace','skutečná poloha lodě', 'začátek pozorování', 'Location', 'best');
title('Pozice lodi');
xlabel('x'); ylabel('y');


% reziduum pozice  

figure(3); hold on;
plot(t, kfres.resid(2,t), 'c-', 'Linewidth', 2);
plot(t, kfres.resid(1,t),'b-','Linewidth', 2);
legend('predikované residuum', 'filtrované residuum');
title('Residuum pozice');
xlabel('t')

mse_filt_x  = (sim.x(1, 1:end-1) - kfres.xf(1, :)) * (sim.x(1, 1:end-1) - kfres.xf(1, :))';
mse_filt_y  = (sim.x(2, 1:end-1) - kfres.xf(2, :)) * (sim.x(2, 1:end-1) - kfres.xf(2, :))';
mse = mse_filt_x + mse_filt_y 

% reziduum rychlosti
cartesian = zeros(2, N);
velocity_car = zeros(1, N);

for ind = 2:N
    cartesian(1,ind) = atand(z(1, ind)/z(1,ind-1))*cos(sqrt(z(1,ind-1)^2 +z(1,ind)^2));
    cartesian(2,ind) = atand(z(1, ind)/z(1,ind-1))*sin(sqrt(z(1,ind-1)^2 +z(1,ind)^2));
    velocity_car(1,ind) = sqrt((cartesian(1,ind)-cartesian(1,ind-1))^2 + (cartesian(2,ind)-cartesian(2,ind-1))^2);

    cartesian_pred(1,ind) = atand(kfres.xp(1, ind)/kfres.xp(1,ind-1))*cos(sqrt(kfres.xp(1,ind-1)^2 +kfres.xp(1,ind)^2));
    cartesian_pred(2,ind) = atand(kfres.xp(1, ind)/kfres.xp(1,ind-1))*sin(sqrt(kfres.xp(1,ind-1)^2 +kfres.xp(1,ind)^2));
    velocity_car_pred(1,ind) = sqrt((cartesian_pred(1,ind)-cartesian_pred(1,ind-1))^2 + (cartesian_pred(2,ind)-cartesian_pred(2,ind-1))^2);
end 


figure(4)
plot(t, velocity_car(1,t), 'r-')
hold on 
plot(t, velocity_car_pred(1,t), 'b-')

figure(5);
diff = abs(velocity_car_pred -velocity_car);
plot(t, diff(1,t))

% stav: [poloha, rychlost]  
% dynamics simulation
% function [x_new, u_new] = state_eq(x, u, w, wu)
%     global f B
%     % u neznamo, v simulaci nahodna prochazka 
%     u_new = eye(2)*u + wu;
%     x_new = f*x + B*u_new + w;
% end
% % measurement simulation --- get new observation
% function new_z = get_observ(x, v)
%    Cx = [atan2(x(2), x(1)); sqrt(x(2)^2 + x(1)^2)];
%    new_z = Cx + v;
% end


function h = getH(x)
    h = [
            atan2(x(2), x(1));
            sqrt(x(2)^2 + x(1)^2);
        ];
end
% jacobian h based on x1, x2
function H = getHJacob(x) 
    H = zeros(2,6); %dim (R,Q)
    H(1:2,1:2) = [
            -x(2)/(x(1)^2+x(2)^2),       x(1)/(x(1)^2 + x(2)^2);
             x(1)/sqrt(x(1)^2 + x(2)^2), x(2)/sqrt(x(1)^2+x(2)^2)
        ];
end


% KALMAN filtr
function res = kf(N, z)
    global P0 m_x0
    % init result struct
    res.xp(:,1) = mvnrnd(m_x0, P0, 1)'; % init x0 -- nejistota poc. odhadu
    res.pp{1} = P0;
    for k = 1:1:N-1
        [res.xf(:,k), res.pf{k}, res.resid(:,k)] = kff(res.xp(:,k), res.pp{k}, z(:,k));
        [res.xp(:, k+1), res.pp{k+1}] = kfp(res.xf(:,k), res.pf{k});
    end
end
function [xp, pp] = kfp(xf, pf)
    global F Q
    xp = F*xf;
    pp = F*pf*F' + Q;
end
function [xf, pf, resid] = kff(xp, pp, z)
    global R
    H = getHJacob(xp);
    h = getH(xp);
    K  = pp * H' * (H*pp*H' + R)^(-1);
    xf = xp + K * (z-h); % bacha v reziduu h() ne jacob. H
    pf = pp - K *H*pp;
    resid = z-h;
end



% % multivariate gaussian sampling 
% function y = mvgrnd(m,sigma,n)
%     % m is the mean vector, sigma is the variance-covariance matrix and n is the number of iterations
% 
%     U = chol(sigma); 
%     % Cholesky decomposition. U is an upper triangular matrix
% 
%     d = length(m);
% 
%     for i = 1:n
%         y(i,1:d) = m' + randn(1,d)*U;
%     end
% end