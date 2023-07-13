function [theta, P, epsilon] = RPEM(fi, P_prev, theta_prev, y, lambda)
    yF = filter([1 0], [1 theta_prev(3)], fi(1)); %filtrovany y(t-1)
    uF = filter([1 0], [1 theta_prev(3)], fi(2)); %filtrovany u(t-1)
    epsilonF = filter([1 0], [1 theta_prev(3)], fi(3)); %filtrovany epsilon(t-1)
    psi = [yF uF epsilonF]'; %novy vektor regeresoru   
    epsilon = y - fi'*theta_prev; %chyba odhadu
    L = (P_prev * psi)/(lambda + psi'*P_prev*psi); %zisk     
    P = (P_prev - (P_prev*psi*psi'*P_prev)/(lambda + psi'*P_prev*psi))/lambda; %nova kovariancni matice
    theta = theta_prev + L*epsilon; %novy odhad
end