% rekurzivni MNC a IVM

clc; clear all; close all


load("data.mat")

y = data(:,1);
u = data(:,2);
T = size(y,1);

zobr1 = 300;
zobr2 = 1000000;

% nejmenší čtverce 
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

figure(1); 
plot(1:zobr1, -theta(1,1:zobr1)',1:zobr1, theta(2,1:zobr1)','LineWidth',2);
hold on 
plot(1:zobr1, ones(1,zobr1)*0.6 ,"k",1:zobr1, ones(1,zobr1)*1,"k");
legend('a', 'b')
title("rekurzivni MNC")
xlabel("iterace")
ylabel("odhad parametrů")

figure(2); 
plot([1,2, 3, 4, 5], [-theta(1,10^2),-theta(1,10^3),-theta(1,10^4),-theta(1,10^5),-theta(1,10^6)],"*",...
    [1,2, 3, 4, 5], [theta(2,10^2),theta(2,10^3),theta(2,10^4),theta(2,10^5),theta(2,10^6)],"*",'LineWidth',2);
hold on 
plot(1:5, ones(1,5)*0.6 ,"k",1:5, ones(1,5)*1,"k");
legend('a', 'b')
title("rekurzivni MNC")
xlabel("iterace")
ylabel("odhad parametrů")
xticks([1,2, 3, 4, 5])
xticklabels({'10^2', '10^3', '10^4', '10^5', '10^6'})
ylim([0.5 1.2])

%% rekurzivni prid. promena -- dodani xi, S2 IVM

P = eye(2) * 1000000; % init P - duvera apriorni informaci

theta1 = zeros(2,T);
theta1(:, 2) = [-2,2]; %apriorni info
for t = 3:T
   xi = [u(t-1); u(t-2)];
   phi = [-y(t-1); u(t-1)]; 

   e = y(t) - phi' * theta1(:, t-1) ; %chyba predikce -- inovace
   L = (P * xi) / (phi' * P * xi + 1);
   theta1(:,t) = theta1(:, t-1) + L*e;
   P = (eye(2) - L*phi')*P; %update P 
end

theta1(:, end)

figure(3); 
plot(1:zobr1, -theta1(1,1:zobr1)',1:zobr1, theta1(2,1:zobr1)','LineWidth',2);
hold on 
plot(1:zobr1, ones(1,zobr1)*0.6 ,"k",1:zobr1, ones(1,zobr1)*1,"k");
legend('a', 'b')
title("rekurzivni IVM")
xlabel("iterace")
ylabel("odhad parametrů")
ylim([0 2])

figure(4); 
plot([1,2, 3, 4, 5], -[theta1(1,10^2),theta1(1,10^3),theta1(1,10^4),theta1(1,10^5),theta1(1,10^6),],"*",...
    [1,2, 3, 4, 5], [theta1(2,10^2),theta1(2,10^3),theta1(2,10^4),theta1(2,10^5),theta1(2,10^6),],"*",'LineWidth',2);
hold on 
plot(1:5, ones(1,5)*0.6 ,"k",1:5, ones(1,5)*1,"k");
legend('a', 'b')
title("rekurzivni IVM")
xlabel("iterace")
ylabel("odhad parametrů")
xticks([1,2, 3, 4, 5])
xticklabels({'10^2', '10^3', '10^4', '10^5', '10^6'})
ylim([0.5 1.2])


%% pseudometoda rozsirene MNC -> odhad i c

P = eye(3) * 1000000; % init P - duvera apriorni informaci

theta2 = zeros(3,T);
e = zeros(T,1);
theta2(:, 1) = [-2,2,2]; %apriorni info -- pridani c .. sumu
% C nevolit nulove?!

for t = 2:T
   phi = [-y(t-1); u(t-1); e(t-1)]; 
   e(t) = y(t) - phi' * theta2(:, t-1); %chyba predikce -- inovace
   L = (P * phi) / (phi' * P * phi + 1);
   theta2(:,t) = theta2(:, t-1) + L*e(t);
   P = (eye(3) - L*phi')*P; %update P 
end

theta2(:, end)

figure(5); 
plot(1:zobr1, -theta2(1,1:zobr1)',1:zobr1, theta2(2,1:zobr1)',1:zobr1, theta2(3,1:zobr1),'LineWidth',2);
hold on 
plot(1:zobr1, ones(1,zobr1)*0.6 ,"k",1:zobr1, ones(1,zobr1)*1,"k",1:zobr1, ones(1,zobr1)*0.5,"k");
legend('a', 'b', 'c')
title("rekurzivní rozšiřené MNČ")
xlabel("počet dat")
ylabel("odhad parametrů")
ylim([0 2])

figure(6); 
plot([1,2, 3, 4, 5], -[theta2(1,10^2),theta2(1,10^3),theta2(1,10^4),theta2(1,10^5),theta2(1,10^6),],"*",...
    [1,2, 3, 4, 5], [theta2(2,10^2),theta2(2,10^3),theta2(2,10^4),theta2(2,10^5),theta2(2,10^6),],"*",...
    [1,2, 3, 4, 5], [theta2(3,10^2),theta2(3,10^3),theta2(3,10^4),theta2(3,10^5),theta2(3,10^6),],"*",'LineWidth',2);
hold on 
plot(1:5, ones(1,5)*0.6 ,"k",1:5, ones(1,5)*1,"k",1:5,ones(1,5)*0.5,"k");
legend('a', 'b','c')
title("rekurzivní rozšiřené MNČ")
xlabel("iterace")
ylabel("odhad parametrů")
xticks([1,2, 3, 4, 5])
xticklabels({'10^2', '10^3', '10^4', '10^5', '10^6'})



%% pseudometoda rozsirene MNC -> odhad i c PEM

th=[-2;2;2];
th_all=th;
P=0.9*eye(3);
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

figure(7); 
plot(1:zobr1, -th_all(1,1:zobr1)',1:zobr1, th_all(2,1:zobr1)',1:zobr1, th_all(3,1:zobr1),'LineWidth',2);
hold on 
plot(1:zobr1, ones(1,zobr1)*0.6 ,"k",1:zobr1, ones(1,zobr1)*1,"k",1:zobr1, ones(1,zobr1)*0.5,"k");
legend('a', 'b', 'c')
title("rekurzivní PEM")
xlabel("počet dat")
ylabel("odhad parametrů")
ylim([0 2])

figure(8); 
plot([1,2, 3, 4, 5], -[th_all(1,10^2),th_all(1,10^3),th_all(1,10^4),th_all(1,10^5),th_all(1,10^6),],"*",...
    [1,2, 3, 4, 5], [th_all(2,10^2),th_all(2,10^3),th_all(2,10^4),th_all(2,10^5),th_all(2,10^6),],"*",...
    [1,2, 3, 4, 5], [th_all(3,10^2),th_all(3,10^3),th_all(3,10^4),th_all(3,10^5),th_all(3,10^6),],"*",'LineWidth',2);
hold on 
plot(1:5, ones(1,5)*0.6 ,"k",1:5, ones(1,5)*1,"k",1:5,ones(1,5)*0.5,"k");
legend('a', 'b','c')
title("rekurzivní PEM")
xlabel("iterace")
ylabel("odhad parametrů")
xticks([1,2, 3, 4, 5])
xticklabels({'10^2', '10^3', '10^4', '10^5', '10^6'})


