function [theta, P, epsilon] = RMSE_extended(fi, P_prev, theta_prev, y)
    epsilon = y - fi'*theta_prev; %chyba odhadu     
    L = (P_prev * fi)/(1 + fi'*P_prev*fi); %zisk
    P = P_prev - ((P_prev*fi*fi'*P_prev)/(1 + fi'*P_prev*fi)); %nova kovariancni matice
    theta = theta_prev + L*epsilon; %novy odhad
end