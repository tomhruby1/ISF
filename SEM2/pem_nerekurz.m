clc; clear all;close all;

load("data.mat")

y = data(:,1);
u = data(:,2);

N = [10^2, 10^3, 10^4, 10^5, 10^6];
pem = zeros(3,5);


for cycle = 1:5
    cycle

theta=[1;1;1];
thetaN=zeros(3,1);


while norm(theta-thetaN)>1e-5
    
    thetaN=theta;
    sum1=zeros(3,3);
    sum2=zeros(3,1);
    eps=0;
    psi=[0 0 0]';
    for i=2:N(cycle)
        sum1=sum1+psi*psi';
        sum2=sum2+psi*eps;
        psi=-thetaN(3)*psi+[-y(i-1) u(i-1) eps]';
        eps=y(i)-[-y(i-1) u(i-1) eps]*thetaN;

    end
    theta=thetaN+inv(sum1)*sum2; 

end
pem(:,cycle) = theta
end

pem
