clc;
clear all;

%% 3., 4. pouziti rekurzivnich metod
N = 1e6; %max. pocet kroku

%Parametry systemu
a = 0.6;
b = 1;
c = 0.5;
lambda = 1; %variance sumu

%Promenne systemu
u = randn(N, 1); %vstupni signal - bily sum
y = zeros(N, 1); %vystup systemu
y(1, :) = 0; %vystup systemu - pocatecni podminka
e = randn(N, 1)*lambda; %sumy

%Promenne odhadu parametru systemu pomoci rozsirene metody nejmensich ctvercu
theta_MSE = zeros(3, N); %odhadovane parametry
theta_MSE(:, 1) = [0.1 0.1 0.1]'; %odhadovane parametry - pocatecni hodnota
P_MSE = zeros(3, 3, N); %kovarinacni matice chyby odhadu                                                            
P_MSE(:, :, 1) = diag([1 1 1]); %kovarinacni matice chyby odhadu - pocatecni hodnota
epsilon_MSE = zeros(1, N); %chyba odhadu
epsilon_MSE(:, 1) = 0; %chyba odhadu - pocatecni hodnota
          
%Promenne odhadu parametru systemu pomoci rekurzivni metody pridavne promenne
theta_IVM = zeros(2, N); %odhadovane parametry
theta_IVM(:, 2) = [0.1 0.1]'; %odhadovane parametry - pocatecni hodnota
P_IVM = zeros(2, 2, N); %kovarinacni matice chyby odhadu                                                            
P_IVM(:, :, 2) = diag([1 1]); %kovarinacni matice chyby odhadu - pocatecni hodnota
              
%Promenne odhadu parametru systemu pomoci rekurzivni metody chyby predikce
theta_PEM = zeros(3, N); %odhadovane parametry
theta_PEM(:, 1) = [0.1 0.1 0.1]';%odhadovane parametry - pocatecni hodnota
P_PEM = zeros(3, 3, N); %kovarinacni matice chyby odhadu                                                            
P_PEM(:, :, 1) = diag([1 1 1]); %kovarinacni matice chyby odhadu - pocatecni hodnota
lambda_PEM = 0.995; %koeficient exponencialniho zapominani dat
epsilon_PEM = zeros(1, N); %chyba odhadu
epsilon_PEM(:, 1) = 0; %chyba odhadu - pocatecni hodnota

%Generovani dat a rekurzivni odhady paremetru systemu
for t = 2:N
    y(t) = a*y(t-1) + b*u(t-1) + e(t) + c*e(t-1); %krok systemu   
    [theta_MSE(:, t), P_MSE(:, :, t), epsilon_MSE(:, t)] = RMSE_extended([y(t-1) u(t-1) epsilon_MSE(t-1)]', P_MSE(:, :, t-1), theta_MSE(:, t-1), y(t)); %rekurzivni metoda nejmensich ctvercu
    [theta_PEM(:, t), P_PEM(:, :, t), epsilon_PEM(:, t)] = RPEM([y(t-1) u(t-1) epsilon_PEM(t-1)]', P_PEM(:, :, t-1), theta_PEM(:, t-1), y(t), lambda_PEM); %rekurzivni metoda chyby predikce
    if t > 2
        [theta_IVM(:, t), P_IVM(:, :, t)] = RIVM([y(t-1) u(t-1)]', [u(t-2) u(t-1)]', P_IVM(:, :, t-1), theta_IVM(:, t-1), y(t)); %rekurzivni metoda pridavne promenne
    end
end

%Porovnani vysledku s realnymi hodnotami parametru
theta_real = [a b c]'
theta_RMSE_extended = theta_MSE(:, N)
theta_RIVM = theta_IVM(:, N)
theta_RPEM = theta_PEM(:, N)

%% 3., 4. vykresleni grafu
k = 1:1:N;
real_par = [a b c]'*ones(1, N);
 
%% 3., 4. rozsirena rekurzivni metoda nejmensich ctvercu
for p = 1:3
    figure;
    hold on;
    grid on;
    ax = gca;
    ax.FontSize = 14;
    plot(k, real_par(p, :), 'b');
    plot(k, theta_MSE(p, :), 'r');
    legend('Skuteèný parametr', 'Odhad')
    xlabel('Èas');
    ylabel('Hodnota parametru');
    
    switch p
        case 1
            title('Rekurzivní odhad a rozšíøenou metodou nejmenších ètvercù');
             saveas(gcf,'RMSE_a','epsc');
        case 2
            title('Rekurzivní odhad b rozšíøenou metodou nejmenších ètvercù');
             saveas(gcf,'RMSE_b','epsc');
        case 3
            title('Rekurzivní odhad c rozšíøenou metodou nejmenších ètvercù');
             saveas(gcf,'RMSE_c','epsc');
    end  
end

%% 3., 4. rekurzivni metoda pridavne promenne
for p = 1:2
    figure;
    hold on;
    grid on;
    ax = gca;
    ax.FontSize = 14;
    plot(k, real_par(p, :), 'b');
    plot(k, theta_IVM(p, :), 'r');
    legend('Skuteèný parametr', 'Odhad')
    xlabel('Èas');
    ylabel('Hodnota parametru');
    
    switch p
        case 1
            title('Rekurzivní odhad a metodou pøídavné promìnné');
            saveas(gcf,'RIVM_a','epsc');
        case 2
        	title('Rekurzivní odhad b metodou pøídavné promìnné');
            saveas(gcf,'RIVM_b','epsc');
    end  
end

%% 3., 4. rekurzivni metoda chyby predikce
for p = 1:3
    figure;
    hold on;
    grid on;
    ax = gca;
    ax.FontSize = 14;
    plot(k, real_par(p, :), 'b');
    plot(k, theta_PEM(p, :), 'r');
    legend('Skuteèný parametr', 'Odhad')
    xlabel('Èas');
    ylabel('Hodnota parametru');
    
    switch p
        case 1
            title('Rekurzivní odhad a metodou chyby predikce');
            saveas(gcf,'RPEM_a','epsc');
        case 2
        	title('Rekurzivní odhad b metodou chyby predikce');
            saveas(gcf,'RPEM_b','epsc');
        case 3
            title('Rekurzivní odhad c metodou chyby predikce');
            saveas(gcf,'RPEM_c','epsc');
    end  
end

%% 3., 4. porovnani rekurzivnich odhadu
for p = 1:3
    figure;
    hold on;
    grid on;
    ax = gca;
    ax.FontSize = 14;    
    plot(k, real_par(p, :), 'k--');   
    if p < 3
        plot(k, theta_IVM(p, :), 'b');
    end
    plot(k, theta_MSE(p, :), 'r');
    plot(k, theta_PEM(p, :), 'g');
    if p < 3
        legend('Reálné parametry', 'IVM', 'MSE', 'PEM');
    else
        legend('Reálné parametry', 'MSE', 'PEM');
    end
    xlabel('Èas');
    ylabel('Hodnota parametru');
    
    switch p
        case 1
            title('Porovnaní rekurzivních odhadù a');
             saveas(gcf,'Rcomp_a','epsc');
        case 2
        	title('Porovnaní rekurzivních odhadù b');
             saveas(gcf,'Rcomp_b','epsc');
        case 3
            title('Porovnaní rekurzivních odhadù c');
             saveas(gcf,'Rcomp_c','epsc');
    end  
end

%% 3., 4. nerozumna pocatecni podminka a porovnani rekurzivniho rozsireneho MSE a PEM

%Promenne odhadu parametru systemu pomoci rozsirene metody nejmensich ctvercu
theta_MSE_WIC = zeros(3, N); %odhadovane parametry
theta_MSE_WIC(:, 1) = [100 100 100]'; %odhadovane parametry - pocatecni hodnota
P_MSE_WIC = zeros(3, 3, N); %kovarinacni matice chyby odhadu                                                            
P_MSE_WIC(:, :, 1) = diag([1 1 1]); %kovarinacni matice chyby odhadu - pocatecni hodnota
epsilon_MSE_WIC = zeros(1, N); %chyba odhadu
epsilon_MSE_WIC(:, 1) = 0; %chyba odhadu - pocatecni hodnota        
              
%Promenne odhadu parametru systemu pomoci rekurzivni metody chyby predikce
theta_PEM_WIC = zeros(3, N); %odhadovane parametry
theta_PEM_WIC(:, 1) = [100 100 100]';%odhadovane parametry - pocatecni hodnota
P_PEM_WIC = zeros(3, 3, N); %kovarinacni matice chyby odhadu                                                            
P_PEM_WIC(:, :, 1) = diag([1 1 1]); %kovarinacni matice chyby odhadu - pocatecni hodnota
lambda_PEM_WIC = 0.995; %koeficient exponencialniho zapominani dat
epsilon_PEM_WIC = zeros(1, N); %chyba odhadu
epsilon_PEM_WIC(:, 1) = 0; %chyba odhadu - pocatecni hodnota


%Generovani dat a rekurzivni odhady paremetru systemu pri spatne pocatecni podmince
for t = 2:N
    y(t) = a*y(t-1) + b*u(t-1) + e(t) + c*e(t-1); %krok systemu   
    [theta_MSE_WIC(:, t), P_MSE_WIC(:, :, t), epsilon_MSE_WIC(:, t)] = RMSE_extended([y(t-1) u(t-1) epsilon_MSE_WIC(t-1)]', P_MSE_WIC(:, :, t-1), theta_MSE_WIC(:, t-1), y(t)); %rekurzivni metoda nejmensich ctvercu
    [theta_PEM_WIC(:, t), P_PEM_WIC(:, :, t), epsilon_PEM_WIC(:, t)] = RPEM([y(t-1) u(t-1) epsilon_PEM_WIC(t-1)]', P_PEM_WIC(:, :, t-1), theta_PEM_WIC(:, t-1), y(t), lambda_PEM_WIC); %rekurzivni metoda chyby predikce   
end

%Porovnani vysledku s realnymi hodnotami parametru
theta_real = [a b c]'
theta_RMSE_extended_WIC = theta_MSE_WIC(:, N)
theta_RPEM_WIC = theta_PEM_WIC(:, N)

%% 3., 4. porovnani rekurzivnich odhadu pri spatne pocatecni podmince
for p = 1:3
    figure;
    hold on;
    grid on;
    ax = gca;
    ax.FontSize = 14;   
    plot(k, real_par(p, :), 'k--');
    plot(k, theta_MSE_WIC(p, :), 'r');
    plot(k, theta_PEM_WIC(p, :), 'g');
    legend('Reálné parametry', 'MSE', 'PEM');
    xlabel('Èas');
    ylabel('Hodnota parametru');
    
    switch p
        case 1
            title('Porovnaní rekurzivních odhadù a');
             saveas(gcf,'Rcomp_a_WIC','epsc');
        case 2
        	title('Porovnaní rekurzivních odhadù b');
             saveas(gcf,'Rcomp_b_WIC','epsc');
        case 3
            title('Porovnaní rekurzivních odhadù c');
             saveas(gcf,'Rcomp_c_WIC','epsc');
    end  
end
