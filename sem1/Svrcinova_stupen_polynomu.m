clear all
clc
close all

%na�ten� dat

t = xlsread("Svrcinova_barc_rocni.xlsx","A2:A104");
y = xlsread("Svrcinova_barc_rocni.xlsx", "B3:B105");

%pomocn� v�po�ty a vykreslen� dat
n=length(y); %po�et pozorov�n�
maxp = 5; %maxim�ln� stupe� polynomu

%%

%odhad stupn� regresn�ho polynomu pomoci metody nejmen��ch �tverc�
%matlab m� pro toto fuknci
%yodh = polyval(polyfit(t,y,p));
%my to ale naprogramujeme:

for p = 1:maxp %projede v�echny stupn� polynomu
    
    X = ones(n,1); %vytvo��me matici X
    for i = 1:p
        X =[X, t.^i];
    end
 
    betaodhad = X'*X\(X'*y); %odhad parametr� metodou nejmen��ch �tverc�
    yodhad = X*betaodhad; %odhad vyrovn�van�ch hodnot
    krit(p) = norm(y-yodhad)^2/(n-p); %kritick� hodnota
    
    
    %kriteria pro odhad stupn� polynomu
    s2(p) = krit(p);
    AIC(p) = log(s2(p))+2*p/n;
    BIC(p) = n* log(s2(p)) + p*log(n);
    FPE(p) = s2(p)*(1+2*p/(n-p));
    R2adj(p)= 1-s2(p)/var(y);
    R2(p) = 1-s2(p)*(n-p)/var(y)*(n-1);
    
end

%najdu nejlep�� hodnotu kriteria a vyp�u stupen polynomu
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

%v�po�et a zabrazen� kritick� hodnoty
% plot(krit,'o')
% grid on
% ylim([0,1])

%vykreslen� kriterii
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


%vykreslen� p�vodn�ch dat
figure
plot(t,y,":o") %vykreslen� dat
grid on %m��ka
hold on %propoj� s dal��mi obr�zky

%odhad parametr� a vyrovnan�ch hodnot pro opti�ln� stupe� polynomu
[b,S] = polyfit(t,y,p_opt);
y_odhad2 = polyval(b,t);
%vykreslen� polynomu
plot(t,y_odhad2,'.-','color','g')
hold on

%predikce na p��t�ch 5 let
h = 5;
plot(max(t)+1:max(t) + h,polyval(b,max(t)+1:max(t) + h),'x-','color','r')
xlabel('Cas t')
ylabel('Prumerne rocni teploty')
title(['Nejvhodn�j�� stupe� polynomu je ',num2str(p_opt)])
hold on

%zobrazen� konfiden�n�ch interval�
[y_odhad3,delta] = polyconf(b,t,S,'simopt','on','predopt','observation');
plot(t,y_odhad3+delta,'b--')
plot(t,y_odhad3-delta,'b--')



