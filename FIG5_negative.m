clear all
clc

format long
global e m a h K 
%alpha1_c=0.0650;%alpha1_c=-0.0868;



Tr=290.15;
e=5.99189*10^(-4);
a=217.59894;h=0.004;m=0.12;K=9;Sr=6;ropt=78.8399;

eps=0.001;
del=0.1;


Tstart=273;Tend=290;

syms R C T al

f=(ropt*exp(-(T-Tr)^2/(2*Sr*Sr)))*R*(1-(R/K))-((a*R*C)/(1+a*h*R));
g=(e*a*R*C)/(1+a*h*R)-m*C;
g2=al;

Cst=((ropt*exp(-(T-Tr)^2/(2*Sr*Sr)))/a)*(1-(R/K))*(1+a*h*R);


s1=diff(f,C);
s2=diff(f,T);
ff1=(g*s1+del*(g2*s2));
f1=subs(ff1,C,Cst);
ff2=-al*del*diff(f,R);
f2=subs(ff2,C,Cst);


al1=-0.1;
f11=subs(f1,al,al1);
f22=subs(f2,al,al1);

eqns=[f11==0,f22==0];
s=solve(eqns,[R,T]);

J=jacobian([f11,f22], [R;T]);

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
xx1=0:0.1:10;
for i=1:length(xx1)
yy1(i)=(m1*(xx1(i)-xo))+yo;
yy2(i)=(m2*(xx1(i)-xo))+yo;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Separting regions with shade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
am=find(yy2>s.T(a1));
am1=am(end);

FP=((a*h*K-1)/(2*a*h));


x1=[FP Rstar xx1(am1+1:end) xx1(end) 10 FP];
y1=[Tend Tstar yy2(am1+1:end) yy2(end) Tend Tend];
fill(x1,y1,'r')
hold on;

%%
x1=[FP xx1(am1+1:end) xx1(end) 10 FP FP];
y1=[Tstar yy2(am1+1:end) yy2(end) Tstart Tstart Tstar];
fill(x1,y1,'y')
hold on;

%%
x1=[0 0 FP FP];
y1=[Tstart Tend Tend Tstart];
fill(x1,y1,'g')
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Maximal canard: Eigenvectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(xx1,yy1,'-b.')
hold on;% % 
plot(xx1,yy2,'-b.')
hold on;% % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Fold line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=linspace(Tstart,Tend,50);

v1=(ropt*exp(-(T-Tr).^2/(2*Sr*Sr)));

R=((a*h*K-1)/(2*a*h)).*ones(1,50);
C=(v1./a).*(1-R/K).*(1+a*h*R);

gg=50;
gg1=1;
plot(R(gg1:gg),T(gg1:gg),'-g.');%fold line
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Folded Singularity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(s.R(1),s.T(1),'ko')%%FS
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rstar=((m)/((a*e-a*m*h)));

ini=[5 273.3];
ts=1:0.01:10000;
al=-0.1;

f=@(t,z)[(z(1)*a*(z(1)*e*ropt*exp(-(z(2) - Tr)^2/(2*Sr^2))*(z(1)/K - 1) - (m*ropt*exp(-(z(2)-Tr)^2/(2*Sr^2))*(z(1)/K - 1)*(z(1)*a*h + 1))/a))/(z(1)*a*h + 1) + (z(1)*al*del*ropt*exp(-(z(2)-Tr)^2/(2*Sr^2))*(z(1)/K-1)*(2*z(2)-2*Tr))/(2*Sr^2);-al*del*((z(1)*ropt*exp(-(z(2)-Tr)^2/(2*Sr^2)))/K +(z(1)*a*h*ropt*exp(-(z(2)-Tr)^2/(2*Sr^2))*(z(1)/K-1))/(z(1)*a*h + 1))]
%f=@(t,z)[(z(1)*a*(z(1)*e*ropt*exp(-(z(2) - Tr)^2/(2*Sr^2))*(z(1)/K - 1) - (m*ropt*exp(-(z(2) - Tr)^2/(2*Sr^2))*(z(1)/K - 1)*(z(1)*a*h + 1))/a))/(z(1)*a*h + 1)+(z(1)*al*del*ropt*exp(-(z(2)-Tr)^2/(2*Sr^2))*(z(1)/K - 1)*(2*z(2) - 2*Tr))/(20*Sr^2);(al*del*((z(1)*ropt*exp(-(z(2) - Tr)^2/(2*Sr^2)))/K + (z(1)*a*h*ropt*exp(-(z(2)- Tr)^2/(2*Sr^2))*(z(1)/K - 1))/(z(1)*a*h + 1)))/10];


[ts,zz]=ode45(f,ts,ini);

v1=(ropt*exp(-(zz(:,2)-Tr).^2/(2*Sr*Sr)));
C=(v1./a).*(1-(zz(:,1)./K)).*(1+a*h.*zz(:,1));

%plot3(zz(1:30000:end-25000,2),zz(1:30000:end-25000,1),C(1:30000:end-25000),'-g.')
hold on;


plot(zz(1:30000:end,1),zz(1:30000:end,2),'-r.','LineWidth',2)
 hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            QSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=linspace(Tstart,Tend,50);
v1=(ropt*exp(-(T-Tr).^2/(2*Sr*Sr)));
Rstar=((m)/((a*e-a*m*h))).*ones(1,50);
plot(Rstar,T,'-r.')

axis([0 10 273 290])