clear all
clc
close all

%načtení dat

load isf_2_data.mat
datax = X1';
datay = Y1';


%%

%vykreslení dat 
        figure 
        plot(datax,datay,":o") %vykreslení dat která jsou k dispozici 
        grid on %mřížka
        hold on 
        xlabel('x')
        ylabel('y')
        title(["Data X1,Y1"])
        hold off 
%%
% rozdělení na identifikační a validační 
datax = datax(10:end,1);
datay = datay(10:end,1);
pocet_dat =length(datax);

% poměr testovacích a trénovacích datinv()
pomer_idval = 0.5;

identifikace = floor(pocet_dat*pomer_idval);
validace = pocet_dat-identifikace;

%sudé hodnoty
t = datax(2:2:end);
y = datay(2:2:end);

%predikuji liché hodnoty
t_pred = datax(1:2:end);
y_pred = datay(1:2:end);

 
%pomocné výpočty a inicializace 
n=length(y); %počet pozorování
maxp = 20; %maximální stupeň polynomu
yodhad = zeros(identifikace, maxp);
ypredikce = zeros(validace, maxp);
chyba = zeros(3, maxp);
%%
msetable = zeros(1,maxp);

for p = 1:maxp %projede všechny stupně polynomu
    p
    
        X = ones(n,1); %vytvoříme matici X
        for i = 1:p
            X =[X, t.^i];
        end
        
        %identifikace 
        betaodhad1 = pinv(X)*y;   %řešení ve smysli nejmenšícch čtverců
        yodhad(1:identifikace,1) = X*betaodhad1;   %odhad vyrovnávaných hodnot

        betaodhad2 = X'*X\(X'*y); %řešení normální soustavy rovnic 
        yodhad(1:identifikace,2) = X*betaodhad2; %odhad vyrovnávaných hodnot
       
        [b,S] = polyfit(t,y,p);
        yodhad(1:identifikace,3) = polyval(b,t);


        Xpr = ones(validace,1); %vytvoříme matici X
        for i = 1:p
            Xpr =[Xpr, t_pred.^i];
        end

        %validace   
        ypredikce(1:validace,1) = Xpr*betaodhad1;
        ypredikce(1:validace,2) = Xpr*betaodhad2;
        ypredikce(1:validace,3)= polyval(b,t_pred);

        %odhad chyby
        chyba(:,p) = [immse(ypredikce(:,1),y_pred) ;immse(ypredikce(:,2),y_pred) ;...
            immse(ypredikce(:,3),y_pred)];


        %vykreslení
        figure 
        plot(datax,datay,":o") %vykreslení dat která jsou k dispozici 
        grid on %mřížka
        hold on 
        plot(t,yodhad(1:identifikace,3),'.-','color','g') %vykreslení polynomu pro identifikaci 
        hold on 
        plot(t_pred,ypredikce(1:validace,3),'x-','color','r')
        xlabel('x')
        ylabel('y')
        title(['Polynom stupně ',num2str(p)])
        hold off     

        msetable(1,p) = mse(ypredikce(:,1),y_pred);
end

msetable