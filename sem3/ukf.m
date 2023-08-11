%% UNCENTED KALMAN FILTER
% inicializace
close all; clear; clc;

dt = 0.1;
N = 500;
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
sim.w = mvgrnd(zeros(size(Q,1),1), Q, N)';
sim.w_u = mvgrnd([0;0], eye(2)*ALPHA_Q, N)'; 
sim.v = mvgrnd([0;0], R, N)'; 

% RANDOM WALK
% sim.x = m_x0(1:4);
% sim.u = [0; 0];
% sim.z = [0; 0];
% for k = 1:1:N-1
%     [sim.x(:,k+1), sim.u(:,k+1)] = state_eq(sim.x(:,k), sim.u(:,k), sim.w(1:4,k), sim.w_u(:,k)); % dyn. sim.
%      sim.z(:,k+1) = get_observ(sim.x(:,k+1), sim.v(:,k+1)); % measurement sim
%end
% figure;
% plot(sim.x')
% legend('x', 'y', 'x*', 'y*');
% title('State');

% POSKYTNUTA DATA
load('isf_4_data.mat');
sim.x = x;
sim.z = z;
sim.u = u;
N = size(sim.x, 2);

%% estimace


tic
 kfres = kf(N, sim.z);
toc


%% RESULTS
k_start = 1; % plot K-first steps
k_end = 120;
t = k_start:k_end;

figure; hold on;
plot(kfres.xp(1,t), kfres.xp(2,t), 'cx');
plot(kfres.xp(1,t), kfres.xp(2,t), 'c-');
plot(kfres.xf(1,t), kfres.xf(2,t), 'bx');
plot(kfres.xf(1,t), kfres.xf(2,t), 'b-');
plot(sim.x(1,t), sim.x(2,t), 'gx');
legend('pos pred','pos pred', 'pos filt', 'pos filt', 'gt - simulation');
title('Position');
xlabel('x'); ylabel('y');

% varF = []
% varF_tot = [,
% varP = []
% varP_tot = []
for k = 1:N-1
    varF(:,k) = diag(kfres.pf{1,k});
    varP(:,k) = diag(kfres.pp{1,k});
    varF_tot(k) = trace(kfres.pf{1,k});
    varP_tot(k) = trace(kfres.pp{1,k});
end

figure; hold on;
plot(t, varP_tot(t));
plot(t, varF_tot(t));
legend( 'tot. pred. var', 'tot. filt. var'); 
title('Variance');

% figure; hold on;
% plot(t, kfres.resid(1,t));
% plot(t, kfres.resid(2,t));
% legend('x resid', 'y resid');
% title('Residuum');


mse_filt_x  = (sim.x(1, 1:end-1) - kfres.xf(1, :)) * (sim.x(1, 1:end-1) - kfres.xf(1, :))';
mse_filt_y  = (sim.x(2, 1:end-1) - kfres.xf(2, :)) * (sim.x(2, 1:end-1) - kfres.xf(2, :))';
mse = mse_filt_x + mse_filt_y 

%%



% stav: [poloha, rychlost]  
% dynamics simulation
function [x_new, u_new] = state_eq(x, u, w, wu)
    global f B
    % u neznamo, v simulaci nahodna prochazka 
    u_new = eye(2)*u + wu;
    x_new = f*x + B*u_new + w;
end
% measurement simulation --- get new observation
function new_z = get_observ(x, v)
   Cx = [atan2(x(2), x(1)); sqrt(x(2)^2 + x(1)^2)];
   new_z = Cx + v;
end


% points sampling
function [chi,w] = uPoints(x,p)
    n = size(x,1);
    % sampling 2n + 1 bodu
    % kappa = 0 (...param rozptylu)
    chi = zeros(2*n+1, n);
    w   = zeros(2*n+1, 1);
    S = chol(p, 'lower'); 
    chi(1, :) = x;
    w(1) = 0;
    for i = 2:n+1
        chi(i, :) = x + sqrt(n) .* S(:,i-1);
    end
    for i = n+2:2*n+1
        chi(i,:) = x - sqrt(n) .* S(:,i-n-1);
    end
    for i = 2:2*n+1
        w(i) = 1/(2*n);
    end
    chi = chi'; % so its [dim x numpoints]
    w = w';
end

function h = getH(x)
    h = [atan2(x(2,:), x(1,:));sqrt(x(1,:).^2 + x(2,:).^2)];
%         h = [
%                 atan2(x(2), x(1));
%                 sqrt(x(2)^2 + x(1)^2);
%             ];
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
    res.xp(:,1) = mvgrnd(m_x0, P0, 1)'; % init x0
    res.pp{1} = P0;
    for k = 1:1:N-1
        [res.xf(:,k), res.pf{k}, res.K{k}] = UKFf(res.xp(:,k), res.pp{k}, z(:,k));
        [res.xp(:, k+1), res.pp{k+1}] = UKFp(res.xf(:,k), res.pf{k});
    end
end
function [xp,Pp]=UKFp(xf,Pf)
    global Q F
    [X,w]=uPoints(xf,Pf);
    Y=F*X;
    xp=Y*w';
    Ydiff=Y-xp;
    Pp=Ydiff.*w*Ydiff'+Q;
end
function [xf,Pf,K]=UKFf(xp,Pp,z)
    global R
    [X,w]=uPoints(xp,Pp);
    Y=getH(X);
    zp=Y*w';
    Ydiff=Y-zp;
    Pxz=(X-xp).*w*Ydiff';
    Pzz=Ydiff.*w*Ydiff'+R;
    dz=z-zp;
    K=Pxz/Pzz;
    xf=xp+K*dz;
    Pf=Pp-K*Pxz';
end




% function [xp, pp] = ukfp(xf, pf)
%     global F Q
%     [X,w] = uPoints(xf, pf);
%     y = F*X; % sampled points through the state eq.
%     xp = y * w';
%     pp = (y-xp) .*w * (y-xp)' + Q;
% end
% function [xf, pf, resid] = ukff(xp, pp, z)
%     global R
%     [X,w] = uPoints(xp, pp);
%     % trasf x do z prostoru ---> z_mean
%     z_transf = zeros(6,13); % 2*
%     for n=1:size(X,2) %over the sample dim
%         % fill in the first 2 dims from h(.)
%         z_transf(1:2,n) = getH(X(:,n)) * w(n);
%     end
%     z_mean = zeros(6,1);
%     z_mean(:) = z_transf * w';
%     pxz = (X - xp) .* w * (X - xp)';
%     R_big = zeros(6,6); R_big(1:2,1:2) = R;
%     pzz = (z_transf - z_mean) .* w * (z_transf - z_mean)' + R_big;
%     z_big = zeros(6,1); z_big(1:2) = z;
%     inov = z_big - z_mean;
%     xf = xp + pxz * pzz^(-1) * inov;
%     pf = pp - pxz * pzz^(-1) * pxz';
% end


% multivariate gaussian sampling 
function y = mvgrnd(m,sigma,n)
    % m is the mean vector, sigma is the variance-covariance matrix and n is the number of iterations

    U = chol(sigma); 
    % Cholesky decomposition. U is an upper triangular matrix

    d = length(m);

    for i = 1:n
        y(i,1:d) = m' + randn(1,d)*U;
    end
end

