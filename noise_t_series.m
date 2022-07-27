%% not working
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

format long

global ropt Tr Sr e a K h sigma m eps
ts=0:0.01:1000; %%time span

sigma=0.08;
Tr=290.15;
e=5.99189*10^(-4);
a=217.59;h=0.004;m=0.12;K=9;Sr=6;ropt=78.8399;
eps=0.01;del=0.1;
T=290;

rr=(ropt*exp(-(T-Tr)^2/(2*Sr*Sr)));

% f=@(t,z)[(1/(eps))*(z(1)*rr*(1-z(1)/K)-(a*z(1)*z(2))/(1+a*h*z(1)));(1/1)*((e*a*z(1)*z(2))/(1+a*h*z(1))-m*z(2))];
% [ts,zz]=ode23t(f,ts,[8 2]); %% ode solver
% 
% %% =============================================
% x0=zz(end,1); 
% y0=zz(end,2);%% setting initial conditions from values obtained using ode solver
Rstar=m/(a*(e-m*h));
Cstar=(rr./a)*(1-Rstar/K)*(1+a*h*Rstar);
x0=Rstar;y0=Cstar;
file_no=0;

mkdir('/home/taran/Desktop/WORK_slow-fast/my work/BINZER PARAM/RATE IN r/noise/DATA/k08')

corr=0;alpha=0.08;
for itr=1%:50
    itr
    file_no=file_no+1;
    fid(file_no) = fopen(sprintf('/home/taran/Desktop/WORK_slow-fast/my work/BINZER PARAM/RATE IN r/noise/DATA/k08/k08_sim%d.dat',file_no),'w' );
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%
    t=0;n=0;
    xd=x0;
    yd=y0;
    dt=0.0001;
    zd=T;
    xi=abs(randn(1,1));
    xid=xi;
    t_max=(313-290)/(alpha);
    n_max=t_max/dt;
    while n<n_max
        n=n+1;
        fprintf(fid(file_no),'%f %f %f %f\n',n,xd,yd,zd);
        
        rr=(ropt*exp(-(zd-Tr)^2/(2*Sr*Sr)));
        
        w=randn(1,1);
        ssd=[(rr*xd*(1-xd/K)) (((a*yd*xd)/(1+xd*a*h))) sigma*sqrt(dt)*w];
        xd=xd+dt*((1/(del*eps))*(ssd(1)-ssd(2)))+ssd(3);
       
        
        yd=yd+dt*((1/del)*((e*a*yd*xd)/(1+xd*a*h)-m*yd));
        zd=zd+dt*alpha;
                
        if xd<0
            xd=0;
        end
        
       
    end
    
end
