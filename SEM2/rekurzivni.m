% rekurzivni MNC a IVM

clc; clear all; 


load("data.mat")

y = data(:,1);
u = data(:,2);
T = size(y,1);


% nejmenší čtvrece 
phi = [-y(1:T-1), u(1:T-1)];
theta_odhad_ls = inv(phi'* phi)*phi' * y(2:T)

%% rekurzivni MNC odhad

P = eye(2) * 10; % init P - duvera apriorni informaci

theta = zeros(2,T);
theta(:, 1) = [-2,2]; %apriorni info
for t = 2:T
   phi = [-y(t-1); u(t-1)]; 
   e = y(t) - phi' * theta(:, t-1); %chyba predikce -- inovace
   L = (P * phi) / (phi' * P * phi + 1);
   theta(:,t) = theta(:, t-1) + L*e;
   P = (eye(2) - L*phi')*P; %update P 
end

theta(:, end)

figure; 
plot(1:T, theta');
legend('a', 'b')
title("rekurzivni MNC")

%% rekurzivni prid. promena -- dodani xi, S2 IVM

P = eye(2) * 1000000; % init P - duvera apriorni informaci

theta = zeros(2,T);
theta(:, 2) = [-2,2]; %apriorni info
for t = 3:T
   xi = [u(t-1); u(t-2)];
   phi = [-y(t-1); u(t-1)]; 
   e = y(t) - phi' * theta(:, t-1); %chyba predikce -- inovace
   L = (P * phi) / (phi' * P * phi + 1);
   theta(:,t) = theta(:, t-1) + L*e;
   P = (eye(2) - L*phi')*P; %update P 
end

theta(:, end)

figure; 
plot(1:T, theta');
legend('a', 'b')
title("rekurzivni IVM")


%% pseudometoda rozsirene MNC -> odhad i c

P = eye(3) * 1000000; % init P - duvera apriorni informaci

theta = zeros(3,T);
e = zeros(T,1);
theta(:, 1) = [-2,2,2]; %apriorni info -- pridani c .. sumu
% C nevolit nulove?!

for t = 2:T
   phi = [-y(t-1); u(t-1); e(t-1)]; 
   e(t) = y(t) - phi' * theta(:, t-1); %chyba predikce -- inovace
   L = (P * phi) / (phi' * P * phi + 1);
   theta(:,t) = theta(:, t-1) + L*e(t);
   P = (eye(3) - L*phi')*P; %update P 
end

theta(:, end)

figure; 
plot(1:T, theta');
legend('a', 'b', 'c')
title("rekurzivní rozšiřené MNČ")

%% pseudometoda rozsirene MNC -> odhad i c PEM

th=[1;1;1];
th_all=th;
P=1*eye(3);
psi=[0 0 0]';
eps=0;
for t=2:length(u)
   psi=-th(3)*psi+[-y(t-1) u(t-1) eps]';
   eps=y(t)-[-y(t-1) u(t-1) eps]*th;
   L=P*psi/(1+psi'*P*psi);
   th=th+L*eps;
   P=(eye(3)-L*psi')*P;
   th_all(:,t)=th;
end

figure; 
plot(1:T, th_all);
legend('a', 'b', 'c')
title("rekurzivní PEM")


