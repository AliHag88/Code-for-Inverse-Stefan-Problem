function [g]=precond_s() % here I don't need to specify input since I am using grad_old everytime
%  t_initial,t_final,Number_of_points, old_gradient
global tmesh

solinit=bvpinit(tmesh,[80 80]);
sol=bvp4c(@bvp4ode,@bvp4bc,solinit);
xint=tmesh;
Sxint=deval(sol,xint);
% plot(xint,Sxint(1,:));
 
g=Sxint(1,:)';

g(1)=0;
 
return     