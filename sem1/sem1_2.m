clc; clear; close all;

load('isf_2_data.mat');

t = X3(1:100)';
y = Y3(1:100)';
X_val = X3(101:end);
Y_val = X3(101:end);

%načtení dat

%t = xlsread("Svrcinova_barc_rocni.xlsx","A2:A104");
%y = xlsread("Svrcinova_barc_rocni.xlsx", "B3:B105");

%pomocné výpočty a vykreslení dat
n=length(y); %počet pozorování
maxp = 5; %maximální stupeň polynomu

%%

%odhad stupně regresního polynomu pomoci metody nejmenších čtverců
%matlab má pro toto fuknci
%yodh = polyval(polyfit(t,y,p));
%my to ale naprogramujeme:

for p = 1:maxp %projede všechny stupně polynomu
    
    X = ones(n,1); %vytvoříme matici X
    for i = 1:p
        X =[X, t.^i];
    end
 
    betaodhad = X'*X\(X'*y); %odhad parametrů metodou nejmenších čtverců
    yodhad = X*betaodhad; %odhad vyrovnávaných hodnot
    krit(p) = norm(y-yodhad)^2/(n-p); %kritická hodnota
    
    
    %kriteria pro odhad stupně polynomu
    s2(p) = krit(p);
    AIC(p) = log(s2(p))+2*p/n;
    BIC(p) = n* log(s2(p)) + p*log(n);
    FPE(p) = s2(p)*(1+2*p/(n-p));
    R2adj(p)= 1-s2(p)/var(y);
    R2(p) = 1-s2(p)*(n-p)/var(y)*(n-1);
    
end

%najdu nejlepší hodnotu kriteria a vypíšu stupen polynomu
%lze urychlit p_opt = find(AIC == min(AIC));

AICMIN = min(AIC);
BICMIN = min(BIC);
FPEMIN = min(FPE);
S2MIN = min(s2);
R2ADJMAX = max(R2adj);
R2MAX = max(R2);

for l = 1:maxp
    if AICMIN == AIC(l)
        disp('AIC min pro polynom stupne')
        p_opt = 1;
        disp(l)
    end
    if BICMIN == BIC(l)
        disp('BIC min pro polynom stupne')
        disp(l)
    end
    if FPEMIN == FPE(l)
        disp('FPE min pro polynom stupne')
        disp(l)
    end
    
end

%výpočet a zabrazení kritické hodnoty
% plot(krit,'o')
% grid on
% ylim([0,1])

%vykreslení kriterii
figure
subplot(2,3,1); plot(1:maxp,AIC,'or');
title('AIC');grid on;
subplot(2,3,2); plot(1:maxp,BIC,'or');
title('BIC');grid on;
subplot(2,3,3); plot(1:maxp,FPE,'or');
title('FPE');grid on;
subplot(2,3,4); plot(1:maxp,s2,'or');
title('s2');grid on;
subplot(2,3,5); plot(1:maxp,R2adj,'or');
title('R2adj');grid on;
subplot(2,3,6); plot(1:maxp,R2,'or');
title('R2');grid on;


%vykreslení původních dat
figure
plot(t,y,":o") %vykreslení dat
grid on %mřížka
hold on %propojí s dalšími obrázky

%odhad parametrů a vyrovnaných hodnot pro optiální stupeň polynomu
[b,S] = polyfit(t,y,p_opt);
y_odhad2 = polyval(b,t);
%vykreslení polynomu
plot(t,y_odhad2,'.-','color','g')
hold on

%predikce na příštích 5 let
h = 5;
plot(max(t)+1:max(t) + h,polyval(b,max(t)+1:max(t) + h),'x-','color','r')
xlabel('Cas t')
ylabel('Prumerne rocni teploty')
title(['Nejvhodnější stupeň polynomu je ',num2str(p_opt)])
hold on

%zobrazení konfidenčních intervalů
[y_odhad3,delta] = polyconf(b,t,S,'simopt','on','predopt','observation');
plot(t,y_odhad3+delta,'b--')
plot(t,y_odhad3-delta,'b--')


