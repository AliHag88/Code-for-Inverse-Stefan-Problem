function dydx=bvp4ode(x,y)

global L_a a_update tmesh



dydx=[y(2) y(1)-interp1(tmesh,a_update, x)/(L_a^2)];

return