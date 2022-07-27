clear all
clc

format long
global e m a h K 
eps=0.001;
del=0.1;
al=0.12;


Tr=290.15;
e=5.99189*10^(-4);
a=217.59894;h=0.004;m=0.12;K=9;Sr=6;ropt=78.8399;


syms R C T 

f=(ropt*exp(-(T-Tr)^2/(2*Sr*Sr)))*R*(1-(R/K))-((a*R*C)/(1+a*h*R));
g=(e*a*R*C)/(1+a*h*R)-m*C;
g2=al;

Cst=((ropt*exp(-(T-Tr)^2/(2*Sr*Sr)))/a)*(1-(R/K))*(1+a*h*R);

s1=diff(f,C);
s2=diff(f,T);
ff1=g*s1+del*(g2*s2);
f1=subs(ff1,C,Cst);
ff2=al*del*diff(f,R);
f2=subs(ff2,C,Cst);

eqns=[f1==0,f2==0];
s=solve(eqns,[R,T]);

J=jacobian([f1,f2], [R;T]);

a1=find(s.T~=0);
Rstar=s.R(a1);
Tstar=s.T(a1)


Jnew=subs(J,{R,T},{Rstar,Tstar});

[eg1 eg2]=eig(Jnew);

xo=Rstar;
yo=Tstar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                    Maximal canard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1=(eg1(2)/eg1(1));
m2=(eg1(4)/eg1(3));
xx1=0:0.1:6;
for i=1:length(xx1)
yy1(i)=(m1*(xx1(i)-xo))+yo;
yy2(i)=(m2*(xx1(i)-xo))+yo;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Maximal canard: Eigenvectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


v1=(ropt*exp(-(yy2-Tr).^2/(2*Sr*Sr)));
C=(v1./a).*(1-(xx1./K)).*(1+a*h.*xx1);

am=find(yy2<s.T(a1));
am1=am(end)+1;
plot3(yy2(am1:5:end),xx1(am1:5:end),C(am1:5:end),'-r.')
%plot3(yy2,xx1,C,'-k.')

hold on; 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ini=[5 302];
ts=1:0.01:10000;
f=@(t,z)[(z(1)*a*(z(1)*e*ropt*exp(-(z(2) - Tr)^2/(2*Sr^2))*(z(1)/K - 1) - (m*ropt*exp(-(z(2)-Tr)^2/(2*Sr^2))*(z(1)/K - 1)*(z(1)*a*h + 1))/a))/(z(1)*a*h + 1) + (z(1)*al*del*ropt*exp(-(z(2)-Tr)^2/(2*Sr^2))*(z(1)/K-1)*(2*z(2)-2*Tr))/(2*Sr^2);-al*del*((z(1)*ropt*exp(-(z(2)-Tr)^2/(2*Sr^2)))/K +(z(1)*a*h*ropt*exp(-(z(2)-Tr)^2/(2*Sr^2))*(z(1)/K-1))/(z(1)*a*h + 1))]
  

[ts,zz]=ode45(f,ts,ini);


v1=(ropt*exp(-(zz(:,2)-Tr).^2/(2*Sr*Sr)));
C=(v1./a).*(1-(zz(:,1)./K)).*(1+a*h.*zz(:,1));

plot3(zz(1:50000:end,2),zz(1:50000:end,1),C(1:50000:end),'-g.')
hold on;