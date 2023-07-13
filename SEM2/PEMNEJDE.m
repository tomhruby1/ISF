clc; clear all; 

% IVM nerekurzivní 

load("data.mat")

y = data(:,1);
u = data(:,2);

%vždy data jen vzít a useknout při novám kodu 
NS = [10^2, 10^3, 10^4] %, 10^5, 10^6];
pocet = size(NS,2); 
ctverce = zeros(2,pocet);
ivm = zeros(2,pocet);
pem = zeros(3,pocet);

for j = 1:pocet
    j
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

    % PEM odhad
    % inicializace paramterů
k = 1000;
theta_hat=[1,1,1];      % počáteční odhad parametrů
alpha = 0.8;



for i = 1:k

    theta_hat0=theta_hat; % na začátku algoritmu předpokládáme, že thety jsou stejné 
    sum1=zeros(3,3);
    sum2=zeros(3,1);

    eps=0;
    psi=[0 0 0]'; %na tyhle deinované nuly se ještšě podíváme
    

    epsilon_filtr = 1;
    y_filtr(1,:) = [0;0;0];

    for i=2:length(y)

        
        y_filtr=-theta_hat(3)* y_filtr + y(i);
        epsilon_filtr= -theta_hat(3)* epsilon_filtr; %tady je ještě nějaké epsilon t ale to já neznám 

        psi= [y_filtr , epsilon_filtr];
        eps= -theta_hat(3)*eps+theta_hat(1)*y(i-1);
    

        sum1=sum1+psi*psi';
        sum2=sum2+eps*psi; %možná vymenit eps a psi?
        
    end

    theta_hat=theta_hat+alpha*inv(sum1)*sum2;
    
end

pem(:, j) = theta_hat;

end