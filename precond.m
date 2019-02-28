function [g]=precond()
%  t_initial,t_final,Number_of_points, old_gradient
global t_ini t_fin Nt 

solinit=bvpinit(linspace(t_ini,t_fin,Nt),[80 80]);
sol=bvp4c(@bvp4ode,@bvp4bc,solinit);
xint=linspace(t_ini,t_fin,Nt);
Sxint=deval(sol,xint);
% plot(xint,Sxint(1,:))
 
g=Sxint(1,:)';
 
return     