clear all
clc
close all

%naètení dat a rozdìlení na identifikaèní a validaèní 

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

%pomocné výpoèty a inicializace 
n=length(y); %poèet pozorování
maxp = 6; %maximální stupeò polynomu
yodhad = zeros(identifikace, maxp);
ypredikce = zeros(validace, maxp);
chyba = zeros(3, maxp);
%%

for p = 1:maxp %projede všechny stupnì polynomu
    
        X = ones(n,1); %vytvoøíme matici X
        for i = 1:p
            X =[X, t.^i];
        end
        
        %identifikace 
        betaodhad1 = pinv(X)*y;   %øešení ve smysli nejmenšícch ètvercù
        yodhad(1:identifikace,1) = X*betaodhad1;   %odhad vyrovnávaných hodnot

        betaodhad2 = X'*X\(X'*y); %øešení normální soustavy rovnic 
        yodhad(1:identifikace,2) = X*betaodhad2; %odhad vyrovnávaných hodnot
       
        [b,S] = polyfit(t,y,p);
        yodhad(1:identifikace,3) = polyval(b,t);


        Xpr = ones(validace,1); %vytvoøíme matici X
        for i = 1:p
            Xpr =[Xpr, t_pred.^i];
        end

        %validace   %DODÌLAT!!!!
        ypredikce(1:validace,1) = Xpr*betaodhad1;
        ypredikce(1:validace,2) = Xpr*betaodhad2;
        ypredikce(1:validace,3)= polyval(b,t_pred);

        %odhad chyby
        chyba(:,p) = [immse(ypredikce(:,1),y_pred) ;immse(ypredikce(:,2),y_pred) ;...
            immse(ypredikce(:,3),y_pred)];


        %vykreslení
        figure 
        plot(datax,datay,":o") %vykreslení dat která jsou k dispozici 
        grid on %møížka
        hold on 
        plot(t,yodhad(1:identifikace,3),'.-','color','g') %vykreslení polynomu pro identifikaci 
        hold on 
        plot(t_pred,ypredikce(1:validace,3),'x-','color','r')
        xlabel('x')
        ylabel('y')
        title(['Polynom stupnì ',num2str(p)])
        hold off         
end

chyba






