function RK02C()

close all
clear all
s = pi/18;
% Solving the ODE for
% 0<t<4 with x1(0)=s and x2(0)=0

%Own RK2 Solver
[t1,x1]=MyRK2Sys(@dxdt,[0 4],[s 0],0.2);                                    %[range of t,x and h]
[t2,x2]=MyRK2Sys(@dxdt,[0 4],[s 0],0.1);
[t3,x3]=MyRK2Sys(@dxdt,[0 4],[s 0],0.05);
%Matlab Solver
[tmat,xmat]=ode23(@dxdt,[0 4],[s 0]);
plot(t1,x1(:,1),'k-',t2,x2(:,1),'b-',t3,x3(:,1),'m-',tmat,xmat(:,1),'ro-')  %[ploting the result]
legend('h=0.2','h=0.1','h=0.05','ode23()');
xlabel('t');
ylabel('x');
title("RK02");

function [t,x]=MyRK2Sys(ODEfunc,tspan,x0,h)
%This function uses 2nd order Runge-Kutta method to solve a
%system of ODEs
t=tspan(1):h:tspan(2);
x(1,:)=x0;
for n=1:length(t)-1
    k1=ODEfunc(t(n),x(n,:))';
    k2=ODEfunc(t(n)+(h/2),x(n,:)+(h/2)*k1)';
    x(n+1,:)=x(n,:)+((1/2)*k1+(1/2)*k2)*h;
end
 
 
function xp=dxdt(t,x)                                                         %[function for the ODE, canonical form]
xp(1)=x(2);
xp(2)=-16.35.*x(1);
xp=xp';