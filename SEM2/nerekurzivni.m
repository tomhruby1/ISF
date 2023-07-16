clc; clear all; 

% IVM nerekurzivní 

load("data.mat")

y = data(:,1);
u = data(:,2);

%vždy data jen vzít a useknout při novám kodu 
NS = [10^2, 10^3, 10^4, 10^5, 10^6];
pocet = size(NS,2); 
ctverce = zeros(2,pocet);
ivm = zeros(2,pocet);
pem = zeros(3,pocet);

for j = 1:pocet
T = NS(j);

    %nejmenší čtverce 
phi = [-y(1:T-1), u(1:T-1)];
theta_odhad_ls = inv(phi'* phi)*phi' * y(2:T);
ctverce(:,j) = theta_odhad_ls;

    % IVM odhad
phi = [-[0; y(1:T-1)], [0;u(1:T-1)]];
epsilon = [[0;u(1:T-1)],[0;0;u(1:T-2)]];
theta_ivm = inv(epsilon' * phi)*epsilon'* y(1:T);
ivm(:, j) = theta_ivm;

end