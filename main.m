
%A sample of how to use Eq(5) in the paper A Dynamics of Surface Temperature Forced by Solar Radiation
%the data used here is from Lubin, D., & Ghiz, M. (2022) "Siple Dome Surface Energy Flux" U.S. Antarctic Program (USAP) Data Center. doi: https://doi.org/10.15784/601540.
%Weiqiang Jing jasonjing@gatech.edu

clear all
load('NSFdata.mat')

global T0
global RnS;
global Q;

global omega0
omega0=2*pi/(24*3600);

%physical parameters for snow
global d; %solar penetration depth
d=4E-3; % d should be ~O(10^-3)
global rho %density
rho=400;
global k %thermal conductivity
k=0.15;
global dt
dt=1;

RnS=@(tau) interp1((1:numel(RnSdata))*1800,RnSdata,tau,'spline','extrap');
Q=@(tau) interp1((1:numel(Qdata))*1800,Qdata,tau,'spline','extrap');
ts=1:dt:numel(Qdata)*1800; %1800s is the sample interval in the original data

%set up ode
%restart every week, update simulation assuming new T0 is observed 


T0=263;
ts1=ts(1:7*24*3600);
[t,Ts1]=ode45(@myfunc,ts1,T0);
%restart simulation assuming new T0 is observed
T0=259;
ts2=ts(7*24*3600+1:14*24*3600);
[t,Ts2]=ode45(@myfunc,ts2,T0);
%the second restart
T0=258;
ts3=ts(14*24*3600+1:end);
[t,Ts3]=ode45(@myfunc,ts3,T0);


Ts=[Ts1;Ts2;Ts3];
Ts=Ts(1:1800:end);

 

figure
% plot(Time,y1)
plot(Time,T_s-273,'blue','LineWidth',1.5);
hold on
 plot(Time,Ts-273,'red--','LineWidth',1.5);
legend('Obs','Dynamic Eqn')
grid on
ylabel('T_s (K)','FontSize', 15,'FontWeight','bold')


figure
syms x
scatter(T_s-273,Ts-273,'b')
R = corr(T_s-273,Ts-273, 'rows','complete');
grid on
hold on
fplot(x,'r','LineWidth',1);
xlim ([-26 -6])
ylim ([-26 -6])
axis square
xlabel('Obs T_s(K)','FontSize', 15,'FontWeight','bold')
ylabel('Modeled T_s(K)','FontSize', 15,'FontWeight','bold')
str=sprintf('r= %1.2f',R);

T = text(min(get(gca, 'xlim'))+4, max(get(gca, 'ylim'))-2, str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');




function dTs=myfunc(t,Ts)
global T0
global omega0
global d;
global k;
global rho
C=2090;
alpha=k/(rho*C);
u=omega0*d^2/alpha;
Is=rho*C*sqrt(alpha);
C1=(1+(u-1)*sqrt(u/2))/(1+u^2);
C2=(u-(u+1)*sqrt(u/2))/(1+u^2);
global dt;
global RnS;
global Q;
dRnSdt=(RnS(t+dt)-RnS(t-dt))/(2*dt);
dTs=sqrt(2*omega0)/Is*(C1*RnS(t)+Q(t)+C2/omega0*dRnSdt)-omega0*(Ts-T0);
end


