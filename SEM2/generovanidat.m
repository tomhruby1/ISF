% generujeme dat 

clear; clc; close all;

a = 0.6;
b = 1; 
c = 0.5;


% počet dat je 501. protože t jde od nuly až do 500
ND = 10^6 + 1; 
t = linspace(0,ND-1,ND);

% definování systému 
y = zeros(ND,1);
y(1,1) = 0;

%chyba 
lambda = 1;
e  =  randn(ND,1)*lambda; % bílý šum

% jako šum je zvolen gauus protož je trvale budící 
u = zeros(ND, 1);
sigma = 1;
u(:,1) = randn(ND, 1)*sigma; % normální

    %generování 
for i = 2:ND
     y(i,1) = a*y(i-1,1) + b*u(i-1,1) + e(i, 1) + c*e(i-1,1);
end

p1 = figure;
scatter(t, y,'.');
title('System');
legend('gaussian');
xlabel("t")
ylabel("y(t)")
%exportgraphics(p1,'system1.pdf')

% data = cat(2,y,u);
% save("data.mat", "data")