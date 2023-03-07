% A
a11 = -0.5; b11 = 0.9; c11 = 0; lambda = 1;
a12 = -0.5; b12 = 0.9; c12 = -0.4;

ND = 501;
t = linspace(0,ND-1,ND);

problem_size = 3;

y1 = zeros(ND,problem_size);
y1(1,1) = 0;
y2 = zeros(ND,problem_size);
y2(1,1) = 0;

sigma = 1;

e  =  randn(ND,1)*lambda; % bílý šum

u = zeros(ND, 3);
for problem = 1:3
    switch problem
        case 1
            u(:,problem)  =  ones(ND,1); % buzení systému
        case 2
            u(:,problem) = zeros(ND,1);
            u(1,problem) = sigma;
        case 3
            u(:,problem) = randn(ND, 1)*sigma;
    end
    
    for i = 2:ND
        y1(i,problem) = -a11*y1(i-1,problem) + b11*u(i-1,problem) + e(i, 1);
        y2(i,problem) = -a12*y2(i-1,problem) + b12*u(i-1,problem) + e(i, 1) + c12*e(i-1,1);
    end
end

figure;
plot(t, y1);
title('System1');
legend('step', 'impulse', 'gaussian');

figure;
plot(t, y2);
title('System2');
legend('step', 'impulse', 'gaussian');

%% B

theta_odhad = zeros(2,3);
for problem = 1:3
    phi = [-y1(1:ND-1,problem), u(1:ND-1,problem)];
    %E = e(1:ND, 1);
    theta_odhad(:,problem) = inv(phi'* phi)*phi' * y1(2:ND,problem); 
end

theta_odhad

%% C??


%%


%theta_odhad = phi' + inv(phi'* phi)*phi'* e(ND-1,1)';
chyba = [a11, b11]'  - theta_odhad;