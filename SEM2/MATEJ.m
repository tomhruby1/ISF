
EPSILON = zeros(1,N);
EPSILON_F = zeros(1,N);
u_F = zeros(1,N);
y_F = zeros(1,N);

EPSILON(1) = y(1); 
EPSILON_F(1) = 1;
u_F(1) = 1;
y_F(1) = 1;
PSI(1,:) = [y_F(1),u_F(1),EPSILON_F(1)];

V_N{k}(1) = EPSILON(1)*EPSILON(1)';
V_N_dot{k}(1) = EPSILON(1)*PSI(1);
V_N_ddot{k}(1) = PSI(1)'*PSI(1);


k=2
while true

    for t=2:N
        EPSILON(t) = y(t) + THETA_hat{k-1}(1)*y(t-1) - THETA_hat{k-1}(2)*u(t-1) - THETA_hat{k-1}(3)*EPSILON(t-1);
        EPSILON_F(t) = EPSILON(t-1) - THETA_hat{k-1}(3)*EPSILON_F(t-1);
        u_F(t) = u(t-1) - THETA_hat{k-1}(3)*u_F(t-1);
        y_F(t) = y(t-1) - THETA_hat{k-1}(3)*y_F(t-1);
        PSI(t,:) = [y_F(t),-u_F(t),-EPSILON_F(t)];
    end
    V_N{k} = EPSILON*EPSILON';
    V_N_dot{k} = EPSILON*PSI;
    V_N_ddot{k} = PSI'*PSI;
    THETA_hat{k}(:) = THETA_hat{k-1}(:) - alpha*inv(V_N_ddot{k})*(V_N_dot{k})';
    
end
