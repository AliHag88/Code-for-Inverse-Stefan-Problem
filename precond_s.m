function [grad]=precond_s(L, s_update, tmesh)

  s_update_fn = @(t)interp1(tmesh, s_update, t);

  % Define ODEs for preconditioning
  ode_s = @(t,y) [y(2), y(1) - s_update_fn(t)/L^2];

  % Define boundary conditions for preconditioning
  bc_s = @(ya, yb) [ya(2) yb(2)];

  % Initialize solver
  solinit = bvpinit(tmesh, [80 80]);

  % Run solver
  sol = bvp4c(ode_s,bc_s,solinit);

  % Evaluate output (preconditioned gradient) on tmesh
  Sint = deval(sol,tmesh);

  % Extract first component and make it a column vector
  grad = Sint(1,:)';

  % Artificially set J_s(0)=0
  grad(1) = 0;
end
