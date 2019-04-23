function dydx=bvp4ode(x,y)

global L_s s_update tmesh



dydx=[y(2) y(1)-interp1(tmesh,s_update, x)/(L_s^2)];

return