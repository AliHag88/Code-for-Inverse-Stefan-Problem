function dydx=bvp4ode(x,y)

global Nt t_ini t_fin L grad

dydx=[y(2) y(1)-interp1(linspace(t_ini,t_fin,Nt),grad, x)/(L^2)];

return






