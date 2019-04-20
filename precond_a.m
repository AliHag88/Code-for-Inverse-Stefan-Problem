function [grad] = precond_a(L, a_update, tmesh)

  a_update_fn = @(t)interp1(tmesh, a_update, t);

  % Define ODEs for preconditioning
  ode_a = @(t,y) [y(2), y(1) - a_update_fn(t)/L^2];

  % Define boundary conditions for preconditioning
  bc_a = @(ya, yb) [ya(2) yb(2)];

  % Initialize solver
  solinit = bvpinit(tmesh, [80 80]);

  % Run solver
  sol = bvp4c(ode_a,bc_a,solinit);

  % Evaluate output (preconditioned gradient) on tmesh
  Sint = deval(sol,tmesh);

  % Extract first component and make it a column vector
  grad = Sint(1,:)';
end
