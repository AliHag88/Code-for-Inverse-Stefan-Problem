clc; clear all; close all;

 len_xmesh=40;
 len_tmesh=30;

 xmesh = linspace(0,1,len_xmesh);  % Space discretization (row)
 tmesh = linspace(0,1,len_tmesh)'; % Time discretization (column)
 
 avals=ones(len_tmesh,1);
 svals=ones(len_tmesh,1);
 
 g=@(t) 2*pi*exp(-4*pi^2*t);
 uInitial=@(x) sin(2*pi*x);

[au_xx_S, u_x_S, u_S, u_T, u]=Forward(xmesh, tmesh, svals, avals, g, uInitial)

u_true=@(x,t) sin(2*pi*x).*exp(-4*pi^2*t);
u__=u_true(xmesh,tmesh);

norm(u-u__);

J=abs(svals(end) - 1)^2 + ...
trapz(xmesh, (u_T - u__(end,:)).^2) + ...
trapz(tmesh, (u_S - u__(:,end)).^2)             
       
c1=1;
B1=1;
k1=1;
tShift=0;

  u_true = @(x,t) c1 + B1*erf(x./(2*sqrt(k1*(t + tShift))));      
       