function new_gradient = precond(tmesh, old_gradient, L)
% precond: Precondition `old_gradient` by projecting it to a Sobolev space.
% Input Arguments:
%    - tmesh: Grid of time values
%    - old_gradient: Previous value of the gradient
%    - L: Preconditioning parameter

%% Initialize solver
yinit = [80, 80];
solinit = bvpinit(tmesh,yinit);

%% Define parameter functions
bvp4bc = @(ya, yb) ...
    deal(...
        ya(2), ...
        yb(2) ...
    );
bvp4ode = @(t, y) ...
    deal(...
        y(2), ...
        y(1) - interp1(tmesh, old_gradient, t) / (L^2) ...
        );

%% Run solver
sol = bvp4c(bvp4ode, bvp4bc, solinit);

%% Extract new gradient
Sxint = deval(sol, xint);
 
new_gradient = Sxint(1,:);
 
end