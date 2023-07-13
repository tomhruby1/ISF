% PEM a rekurzivní PEM - upravit kod 

clear all

load("data.mat")

y = data(:,1);
u = data(:,2);
T = size(y,1);

NS = [10^2, 10^3, 10^4, 10^5, 10^6];
pocet = size(NS,2);

theta=[1;1;1];
thetaN=zeros(3,1);
while norm(theta-thetaN)>1e-5
    thetaN=theta;
    sum1=zeros(3,3);
    sum2=zeros(3,1);
    eps=0;
    psi=[0 0 0]';
    for i=2:length(y)
        sum1=sum1+psi*psi';
        sum2=sum2+psi*eps;
        psi=-thetaN(3)*psi+[-y(i-1) u(i-1) eps]';
        eps=y(i)-[-y(i-1) u(i-1) eps]*thetaN;
    end
    theta=thetaN+inv(sum1)*sum2; 
end

%% 
theta=[1;1;1];
thetaN=zeros(3,1);
while norm(theta-thetaN)>1e-5 
    thetaN=theta;
    sum1=zeros(3,3);
    sum2=zeros(3,1);
    eps=0;
    psi=[0 0 0]';
    for i=2:length(y)
        sum1=sum1+psi*psi';
        sum2=sum2+psi*eps;
        psi=-thetaN(3)*psi+[-y(i-1) u(i-1) eps]';
        eps=y(i)-[-y(i-1) u(i-1) eps]*thetaN;
    end
    theta=thetaN+inv(sum1)*sum2; 
end
theta
%% rekurz

th=[1;1;1];
th_all=th;
P=1*eye(3);
psi=[0 0 0]';
eps=0;
for t=2:length(u)
   psi=-th(3)*psi+[-y(t-1) u(t-1) eps]';
   eps=y(t)-[-y(t-1) u(t-1) eps]*th;
   L=P*psi/(1+psi'*P*psi);
   th=th+L*eps;
   P=(eye(3)-L*psi')*P;
   th_all(:,t)=th;
end
th
plot(1:T,th_all)
