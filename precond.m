function grad = precond(tmesh, L, grad_in)

solinit = bvpinit(tmesh,[80 80]);

ode_s = @(t, y) [y(2) y(1) - interp1(tmesh, grad_in, t)/L^2];

sol = bvp4c(ode_s,@bc_s,solinit);

Sxint = deval(sol,tmesh);

grad = Sxint(1,:)';

grad(1) = 0;
end

function res = bc_s(ya, yb)
    res = [ya(2), yb(2)];
end