function grad = precond(mode, L, tmesh, grad_t, varargin)
  % precond applies "Sobolev preconditioning" to the given arguments
  % - mode is an integer. Options:
  %   mode == 1: Give the output homogeneous Neumann conditions
  %   mode == 2: Give the output a Neumann condition at the RHS and Dirichlet
  %              condition on the LHS.
  % - tmesh is the time grid on which the problem is solved.
  % - L is the weight associated with precondition
  % - grad_t is the original gradient J(t)
  % This function accepts a variable argument list:
  % - If given 5 arguments, we assume that the 5th argument is J(T).

  % Preconditioning doesn't work when the parameter is too small.
  if L < 1e-10
    grad = grad_t;
    return
  end

  % Set gradient at final moment if requested.
  if nargin == 5
    grad_final_moment = varargin{1};
  else
    grad_final_moment = 0;
  end

  % Initialize BVP solver
  solinit = bvpinit(tmesh, [0, 0]);

  % Switch BC setup as set by `mode`
  if mode == 1
    bc_s = @(ya, yb) [ya(2), L*yb(2) - grad_final_moment];
  elseif mode == 2
    bc_s = @(ya, yb) [ya(1), L*yb(2) - grad_final_moment];
  else
    error('Unsupported mode passed to precond script.');
  end

  % Define ODE to be satisfied by the preconditioned gradient
  ode_s = @(t, y) [y(2), (y(1) - interp1(tmesh, grad_t, t))/L];

  % Solve preconditioning problem
  sol = bvp4c(ode_s, bc_s, solinit);

  % Evaluate solution on time grid
  Sxint = deval(sol, tmesh);

  % Extract y_1 (solution) component
  grad = Sxint(1,:)';
end
