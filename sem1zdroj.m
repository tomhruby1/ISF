a = -0.5;
b = 0.9;
c = 0;
 a_odhad = a- c*(1-a^2)/(1+c^2-2*a*c);
 b_odhad = b-b*c*(1-a)/(1+c^2-2*a*c);
 
 
 y_1 = 2.7339;
 y_2 = -0.8919;
 
 a_impuls = a
 b_impuls = (a*y_1+y_2)
 
 %% 
 syms('a', 'b', 'c', 'sigma', 'lambda')
 
m11 = (b^2 * sigma + (1+c^2  - 2*a*c)*lambda)
m12 = 0
m21 = 0 
m22 = sigma

v1 = (-a*b^2*sigma + (c-a) * (1-a*c) * lambda)/(1-a^2)
v2 = b*sigma

M = [m11 m12; m21 m22];
res = inv(M) * [v1; v2]


simplify(res(1))

%%
lambda_real = 1;
sigma_real = 1; 
dosazeno = subs(res, [sigma_real, lambda_real], [1,1])


