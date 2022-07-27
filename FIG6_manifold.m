clc;
clear all ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% surface plot for the manifold, fold line and the QSSA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global e m a h K ropt Tr Sr 
subplot(1,2,1)
Tstart=273;Tend=290; 

Tr=290.15;
e=5.99189*10^(-4);
a=217.59894;h=0.004;m=0.12;K=9;Sr=6;ropt=78.8399;

parts=9;


T=linspace(Tstart,Tend,parts);
R=linspace(0,8,parts);


[T,R]=meshgrid(T,R);
v1=(ropt*exp(-(T-Tr).^2/(2*Sr*Sr)));
C=(v1./a).*(1-(R./K)).*(1+a*h.*R);

parts_2=5;
surf(T(1:parts_2,1:end),R(1:parts_2,1:end),C(1:parts_2,1:end),'FaceAlpha',0.6);
hold on;
surf(T(parts_2:end,1:end),R(parts_2:end,1:end),C(parts_2:end,1:end),'FaceAlpha',0.2);
hold on;

xlabel('T') 
ylabel('R') 
zlabel('C') 
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C1=C(1,:);

yy=C1;


T=linspace(Tstart,Tend,parts);
R=0.*ones(1,parts);



yy2=ones(parts,parts);
for i=1:parts
yy2(i,:)=linspace(0,yy(i),parts);
end

yy3=ones(parts,parts);
for i=1:parts
yy3(i,:)=linspace(yy(i),1,parts);
end
[T,R]=meshgrid(T,R);


surf(T(1:end,1:end),R(1:end,1:end),yy2(1:end,1:end)','FaceAlpha',0.6);
hold on;

surf(T(1:end,1:end),R(1:end,1:end),yy3(1:end,1:end)','FaceAlpha',0.2);
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


T=linspace(Tstart,Tend,parts);


v1=(ropt*exp(-(T-Tr).^2/(2*Sr*Sr)));

R=((a*h*K-1)/(2*a*h)).*ones(1,parts);
C=(v1./a).*(1-R/K).*(1+a*h*R);

 gg=20;
 gg1=1;
 plot3(T(1:end),R(1:end),C(1:end),'-g');%fold line
 hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=linspace(Tstart,Tend,parts);


v1=(ropt*exp(-(T-Tr).^2/(2*Sr*Sr)));

Rstar=((m)/((a*e-a*m*h))).*ones(1,parts);
Cstar=(v1./(a)).*(1-Rstar/K).*(1+a*h*Rstar)
plot3(T,Rstar,Cstar,'-r')
hold on;


 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







subplot(1,2,2)
Tstart=313;Tend=290; 

Tr=290.15;
e=5.99189*10^(-4);
a=217.59894;h=0.004;m=0.12;K=9;Sr=6;ropt=78.8399;

parts=9;


T=linspace(Tstart,Tend,parts);
R=linspace(0,8,parts);


[T,R]=meshgrid(T,R);
v1=(ropt*exp(-(T-Tr).^2/(2*Sr*Sr)));
C=(v1./a).*(1-(R./K)).*(1+a*h.*R);

parts_2=5;
surf(T(1:parts_2,1:end),R(1:parts_2,1:end),C(1:parts_2,1:end),'FaceAlpha',0.6);
hold on;
surf(T(parts_2:end,1:end),R(parts_2:end,1:end),C(parts_2:end,1:end),'FaceAlpha',0.2);
hold on;

xlabel('T') 
ylabel('R') 
zlabel('C') 
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C1=C(1,:);

yy=C1;


T=linspace(Tstart,Tend,parts);
R=0.*ones(1,parts);



yy2=ones(parts,parts);
for i=1:parts
yy2(i,:)=linspace(0,yy(i),parts);
end

yy3=ones(parts,parts);
for i=1:parts
yy3(i,:)=linspace(yy(i),1,parts);
end
[T,R]=meshgrid(T,R);


surf(T(1:end,1:end),R(1:end,1:end),yy2(1:end,1:end)','FaceAlpha',0.6);
hold on;

surf(T(1:end,1:end),R(1:end,1:end),yy3(1:end,1:end)','FaceAlpha',0.2);
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


T=linspace(Tstart,Tend,parts);


v1=(ropt*exp(-(T-Tr).^2/(2*Sr*Sr)));

R=((a*h*K-1)/(2*a*h)).*ones(1,parts)
C=(v1./a).*(1-R/K).*(1+a*h*R)

 gg=20;
 gg1=1;
 plot3(T(1:end),R(1:end),C(1:end),'-g');%fold line
 hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=linspace(Tstart,Tend,parts);


v1=(ropt*exp(-(T-Tr).^2/(2*Sr*Sr)));

Rstar=((m)/((a*e-a*m*h))).*ones(1,parts);
Cstar=(v1./(a)).*(1-Rstar/K).*(1+a*h*Rstar)
plot3(T,Rstar,Cstar,'-r')
hold on;


 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%