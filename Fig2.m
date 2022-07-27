clc;
clear all ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig2 %equilibrium points, manifolds and trajectory evolution for setting T=287
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global e m a h K ropt Tr Sr 



Tr=290.15;
e=5.99189*10^(-4);
a=217.59894;h=0.004;m=0.12;K=9;Sr=6;ropt=78.8399;

T=273;

v1=(ropt*exp(-(T-Tr).^2/(2*Sr*Sr)));

% nullcline
R=linspace(0,K,50);
C=0.*ones(1,50);
plot(R,C,'g.')
hold on;

R=m/(e*a-m*a*h).*ones(1,50);
Cs=(v1./a).*(1-R/K).*(1+a*h*R);
C=linspace(0,1,50);
plot(R,C,'go')%fixed point
%%eqm
R=0;C=0;
plot(R,C,'bo')
hold on;

C=0;R=(K);
plot(R,C,'bo')
hold on;

C=linspace(0,0.4,50);
R=0.*ones(1,50);
plot(R,C,'-r') %manifold

C=linspace(0.4,1,50);
R=0.*ones(1,50);
plot(R,C,'-y') %manifold


for T=287
    T
FP=((a*h*K-1)/(2*a*h));
R=linspace(0,FP,50);

v1=(ropt*exp(-(T-Tr).^2/(2*Sr*Sr)));
C=(v1./a).*(1-(R/K)).*(1+a*h*R);
plot(R,C,'-y') %manifold
hold on;


R=linspace(FP,K,50);

v1=(ropt*exp(-(T-Tr).^2/(2*Sr*Sr)));
C=(v1./a).*(1-(R/K)).*(1+a*h*R);
plot(R,C,'-r') %manifold
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%eqm
R=m/(e*a-m*a*h);
C=(v1./a).*(1-R/K).*(1+a*h*R);
plot(R,C,'bo')
hold on; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%fold point
R=((a*h*K-1)/(2*a*h)).*ones(1,50);
C=(v1./a).*(1-R/K).*(1+a*h*R);
plot(R,C,'-ko')
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%trajectory
ini=[7 0.007];
ts=0:0.01:1000;

f=@(t,z)[(1/eps)*((ropt*exp(-(T-Tr)^2/(2*Sr*Sr))*z(1)*(1-(z(1)/K))-((a*z(1)*z(2))/(1+a*h*z(1)))));((e*a*z(1)*z(2))/(1+a*h*z(1)))-m*z(2)];
[ts,zz]=ode15s(f,ts,ini); %% ode solver

plot(zz(1:end,1),zz(1:end,2),'-k.')
hold on;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
