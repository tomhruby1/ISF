clc;
clear all;
tic

%% 1., 2. nerekurzivni odhady
Ns = [10^2, 10^3, 10^4]

number_of_runs = 10 %kolikrat jedu kaydou metodu pro kazdy N
[s1 s2] = size(Ns)
means_IVM = NaN(s2, 2)
means_PEM = NaN(s2, 3)
vars_IVM = NaN(s2, 2)
vars_PEM = NaN(s2, 3)

for index_N = 1:s2 %pro kazdy N pocet generovanych dat
    
%% 1., 2. parametry systému

    N =  Ns(index_N); %t = 0 ... 501
    a1 = -0.6;
    b1 = 1;
    c1 = 0.5;
    lambda_kva_S = 1;
    a_S = [1, a1];
    b_S = [b1];
    c_S = [1, c1];
    sigma_kva = 1; %variance bileho sumu na vstupu
    u_noise = randn(N, 1)* sqrt(sigma_kva);
    
%% 1., 2. metoda pridavne promenne (nerekurzivni)
    Theta_IVM_results = zeros(2, number_of_runs);


    for i = 1:number_of_runs
    	u = u_noise;
        y = generate_data(a_S, b_S, c_S, u, sqrt(lambda_kva_S), N);
        y_t_minus_1 = [0; y(1:end-1, 1)];
        u_t_minus_1 = [0; u(1:end, 1)];
        u_t_minus_2 = [0;0; u(1:end-1, 1)];
        
        Phi = [-y_t_minus_1 u_t_minus_1];
        Psi = [u_t_minus_1, u_t_minus_2]; %pridavne promenne
        Theta = (inv(Psi' * Phi))*Psi'*y;
        Theta(1) = Theta(1)* -1; %sruktura dle zadani kdzytak odrtranit
        Theta_IVM_results(:,i) = Theta;
    end

    means_IVM(index_N,:) = mean(Theta_IVM_results')
    vars_IVM(index_N,:) =    var(Theta_IVM_results')

    
%% 1., 2. metoda chyby predikce (nerekurzivni)
% hledam optimalni prediktor - jeho parametry jsou parametry systemu
    Theta_PEM_results = zeros(3, number_of_runs); %3 odhaduju i c

    for i = 1:number_of_runs
        u = u_noise;
        y = generate_data(a_S, b_S, c_S, u, sqrt(lambda_kva_S), N);
        y_t_minus_1 = [0; y(1:end-1, 1)];
        u_t_minus_1 = [0; u(1:end, 1)];
        a_ini = 0.2;
        b_ini = 0.2;
        c_ini = 0.2;
        Theta = [a_ini ; b_ini; c_ini];
        alpha = 1;%krok
	
        iteration = 1;
        max_iterations = 50;
        threshold = 0.00001;
    
        while true
            Theta_old = Theta;
            %spocti  epsilon - odhad
            epsilon = filter([1, Theta(1)], [1, Theta(3)],y) - filter([Theta(2), 0], [1, Theta(3)],u_t_minus_1) ;
            epsilon_t_minus_1 = [0; epsilon(1:end-1, 1)];
            %spocti filtrovane yf, uf epsilonf
            yf = filter([1, 0], [1, Theta(3)], y_t_minus_1); 
            uf =  filter([1, 0], [1, Theta(3)], u_t_minus_1);
            epsilonf = filter([1, 0], [1, Theta(3)], epsilon_t_minus_1);

            Psi = +[-yf, +uf, epsilonf];

            subfunction1 =  trace(Psi'*Psi);
            subfunction2 =  Psi'*epsilon;
            % subfunction2
            Theta  = Theta_old + alpha*1/(subfunction1)*(subfunction2);
            criterial_function = norm(Theta - Theta_old);
            if(iteration == max_iterations || criterial_function< threshold)
                break %kontrola konvergence algoritmu/dosayeni max poctu iteraci
            end
                iteration = iteration+1; %pridani iterace    
        end
            Theta(1) = Theta(1)* -1; %sruktura dle zadani kdzytak odrtranit
            Theta_PEM_results(:, i) = Theta;     
    end

    means_PEM(index_N,:) = mean(Theta_PEM_results')
    vars_PEM(index_N,:) =    var(Theta_PEM_results')

end

%% 1., 2. vypis vysledku
means_IVM
means_PEM
vars_IVM
vars_PEM
toc