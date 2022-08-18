clear all
load('NSFdata.mat')


 global T0

T0=261;



global omega0
omega0=2*pi/(24*3600);
global d;
d=4E-3;
global rho
rho=400;
global k
k=0.15;
global dt
dt=1;

global RnS;


RnS=@(tau) interp1((1:numel(RnSdata))*1800,RnSdata,tau,'spline','extrap');

global Q;
Q=@(tau) interp1((1:numel(Qdata))*1800,Qdata,tau,'spline','extrap');
ts=1:dt:numel(Qdata)*1800;







%st up ode

T0=261;
ts1=ts(1:10*24*3600);
[t,Ts1]=ode45(@myfunc,ts1,T_s(1));

T0=258;
ts2=ts(10*24*3600+1:end);
[t,Ts2]=ode45(@myfunc,ts2,T_s(10*48));

Ts=[Ts1;Ts2];

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




return
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


