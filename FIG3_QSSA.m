clear all
clc

format long

global e m a h K ropt Tr Sr 
%% Fig3-QSSA for the slow-fast system%

eps=0.001;
del=0.1;

Tr=290.15;
e=5.99189*10^(-4);
a=217.59894;h=0.004;m=0.12;K=9;Sr=16;ropt=78.8399;


ts=0:0.01:600; %%time span
eps=0.01;

T=273:1:290;
for i=1:length(T)
    i
r=ropt*exp(-(T(i)-Tr)^2/(2*Sr*Sr));
ini=[(m/(e*a-a*h*m)) (r/a)*(1-(m/(e*a-a*h*m))/K)*(1+a*h*(m/(e*a-a*h*m)))];

f=@(t,z)[(1/(eps))*(ropt*exp(-(T(i)-Tr)^2/(2*Sr*Sr))*z(1)*(1-(z(1)/K))-((a*z(1)*z(2))/(1+a*h*z(1))));((e*a*z(1)*z(2))/(1+a*h*z(1)))-m*z(2)];

[ts,zz]=ode23tb(f,ts,ini); %% ode solver
zz1(i)=mean(zz(1:end,1));
zz2(i)=mean(zz(1:end,2));

end

subplot(2,2,1)
zz1=fliplr(zz1);
plot(1:length(zz1),zz1,'-k.')
hold on;
subplot(2,2,2)
zz2=fliplr(zz2);
plot(1:length(zz2),zz2,'-k.')
hold on;

