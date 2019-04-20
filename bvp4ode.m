function dydx=bvp4ode(x,y)

global L s_update tmesh



dydx=[y(2) y(1)-interp1(tmesh,s_update, x)/(L^2)];

return