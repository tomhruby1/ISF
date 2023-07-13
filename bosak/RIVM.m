function [theta, P] = RIVM(fi, ksi, P_prev, theta_prev, y)
    E = y - fi'*theta_prev; %chyba odhadu
    L = (P_prev * ksi)/(1 + fi'*P_prev*ksi); %zisk 
    theta = theta_prev + L*E; %novy odhad
    P = P_prev - ((P_prev*ksi*fi'*P_prev)/(1 + fi'*P_prev*ksi)); %nova kovariancni matice
end