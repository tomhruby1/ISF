clear all
clc
close all

%na�ten� dat a rozd�len� na identifika�n� a valida�n� 

load isf_2_data.mat
datax = X1';
datay = Y1';
pocet_dat =length(datax);
pomer_idval = 0.95;
identifikace = floor(pocet_dat*pomer_idval);
validace = pocet_dat-identifikace;


t = datax(1:identifikace,1);
y = datay(1:identifikace,1);

t_pred = datax(identifikace+1:pocet_dat,1);
y_pred = datay(identifikace+1:pocet_dat,1);

%pomocn� v�po�ty a inicializace 
n=length(y); %po�et pozorov�n�
maxp = 6; %maxim�ln� stupe� polynomu
yodhad = zeros(identifikace, maxp);
ypredikce = zeros(validace, maxp);
chyba = zeros(3, maxp);
%%

for p = 1:maxp %projede v�echny stupn� polynomu
    
        X = ones(n,1); %vytvo��me matici X
        for i = 1:p
            X =[X, t.^i];
        end
        
        %identifikace 
        betaodhad1 = pinv(X)*y;   %�e�en� ve smysli nejmen��cch �tverc�
        yodhad(1:identifikace,1) = X*betaodhad1;   %odhad vyrovn�van�ch hodnot

        betaodhad2 = X'*X\(X'*y); %�e�en� norm�ln� soustavy rovnic 
        yodhad(1:identifikace,2) = X*betaodhad2; %odhad vyrovn�van�ch hodnot
       
        [b,S] = polyfit(t,y,p);
        yodhad(1:identifikace,3) = polyval(b,t);


        Xpr = ones(validace,1); %vytvo��me matici X
        for i = 1:p
            Xpr =[Xpr, t_pred.^i];
        end

        %validace   %DOD�LAT!!!!
        ypredikce(1:validace,1) = Xpr*betaodhad1;
        ypredikce(1:validace,2) = Xpr*betaodhad2;
        ypredikce(1:validace,3)= polyval(b,t_pred);

        %odhad chyby
        chyba(:,p) = [immse(ypredikce(:,1),y_pred) ;immse(ypredikce(:,2),y_pred) ;...
            immse(ypredikce(:,3),y_pred)];


        %vykreslen�
        figure 
        plot(datax,datay,":o") %vykreslen� dat kter� jsou k dispozici 
        grid on %m��ka
        hold on 
        plot(t,yodhad(1:identifikace,3),'.-','color','g') %vykreslen� polynomu pro identifikaci 
        hold on 
        plot(t_pred,ypredikce(1:validace,3),'x-','color','r')
        xlabel('x')
        ylabel('y')
        title(['Polynom stupn� ',num2str(p)])
        hold off         
end

chyba






