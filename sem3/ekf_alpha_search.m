%% EXTENDED KALMAN FILTER
% inicializace
close all; clear; clc;



dt = 0.1;
N = 120;
M = 8;



SEARCH_STEPS = 100;

% define search steps
ALPHA_Q_OG = 1e-3;
ALPHA_P_OG = 1e-6;
aqs = 1e-4:1e-4:0.01; %logspace(1e-4,1e-5, SEARCH_STEPS); %1e-4:1e-4:0.01;
aps = 1e-4:1e-4:0.01; %logspace(1e-7,1e-3, SEARCH_STEPS); %1e-7:1e-4:0.01;
aqs(SEARCH_STEPS+1) = ALPHA_Q_OG;
aqs(SEARCH_STEPS+1) = ALPHA_P_OG;

global f B F Q R m_x0 P0

aq_count = 0;
mse_pred_x = zeros(length(aqs), length(aps));
mse_pred_y = zeros(length(aqs), length(aps));
mse_filt_x = zeros(length(aqs), length(aps));
mse_filt_y = zeros(length(aqs), length(aps));

% GRID SEARCH! 
for ALPHA_Q = aqs
    aq_count = aq_count + 1;
    ap_count = 0;
    for ALPHA_P = aps
        ap_count = ap_count + 1;
        for m = 1:M
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
            Q(5,5) = ALPHA_Q; Q(6,6) = ALPHA_P;

            R = [
                    4e-4 * (pi/180)^2      0;
                    0                    1e-4
                ];

            m_x0 = [20; 50; 0; -12; 0; 0;];


            P0 = eye(6);
            P0(5,5) = ALPHA_P; P0(6,6) = ALPHA_P; %pocatecni podminka pro rizeni


            % precompute noises
            sim.w = mvgrnd(zeros(size(Q,1),1), Q, N)';
            sim.w_u = mvgrnd([0;0], eye(2)*ALPHA_Q, N)'; 
            sim.v = mvgrnd([0;0], R, N)'; 

            sim.x = m_x0(1:4);
            sim.u = [0; 0];
            sim.z = [0; 0];
            for k = 1:1:N-1
                [sim.x(:,k+1), sim.u(:,k+1)] = state_eq(sim.x(:,k), sim.u(:,k), sim.w(1:4,k), sim.w_u(:,k)); % dyn. sim.
                 sim.z(:,k+1) = get_observ(sim.x(:,k+1), sim.v(:,k+1)); % measurement sim
            end

            kfres = kf(N, sim.z);
            
            % calculate MSE and avg over m
%             mse_pred_x(aq_count,ap_count) = mse_pred_x(aq_count,ap_count)+(sim.x(1, :) - kfres.xp(1, :)) * (sim.x(1, :) - kfres.xp(1, :))' /M;
%             mse_pred_y(aq_count,ap_count) = mse_pred_x(aq_count,ap_count)+(sim.x(2, :) - kfres.xp(2, :)) * (sim.x(2, :) - kfres.xp(2, :))' /M;
%             mse_filt_x(aq_count,ap_count) = mse_pred_x(aq_count,ap_count)+(sim.x(1, 1:end-1) - kfres.xf(1, :)) * (sim.x(1, 1:end-1) - kfres.xf(1, :))' /M;
%             mse_filt_y(aq_count,ap_count) = mse_pred_x(aq_count,ap_count)+(sim.x(2, 1:end-1) - kfres.xf(2, :)) * (sim.x(2, 1:end-1) - kfres.xf(2, :))' /M;
            
            % MAE
            mse_pred_x(aq_count,ap_count) = mse_pred_x(aq_count,ap_count)+sum((abs(sim.x(1, :) - kfres.xp(1, :)))) /M;
            mse_pred_y(aq_count,ap_count) = mse_pred_x(aq_count,ap_count)+sum((abs(sim.x(2, :) - kfres.xp(2, :)))) /M;
            mse_filt_x(aq_count,ap_count) = mse_pred_x(aq_count,ap_count)+sum((abs(sim.x(1, 1:end-1) - kfres.xf(1, :)))) /M;
            mse_filt_y(aq_count,ap_count) = mse_pred_x(aq_count,ap_count)+sum((abs(sim.x(2, 1:end-1) - kfres.xf(2, :)))) /M;
        end
        [aq_count, ap_count] 
    end
end
%% 
figure;
the_mse = mse_filt_x + mse_filt_y;

surf(aps, aqs, the_mse)
colorbar
zlim([0, 100000])
zlim


[M,I] = min(the_mse(:));
[I_row, I_col] = ind2sub(size(the_mse),I);
min_idx = [aqs(I_row), aps(I_col)]
M


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
    res.xp(:,1) = mvgrnd(m_x0, P0, 1)'; % init x0 -- nejistota poc. odhadu
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

