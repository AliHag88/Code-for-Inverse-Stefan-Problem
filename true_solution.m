function [u_true_T, u_true_0, u_true_S, aux_true_0, s_star, s_true, u_true, a_true] = true_solution(tmesh, k1, c1, latentHeat, tShift)
  % true_solution: Return parameter functions for "true" solution.
  % See notes.tex for a citation to Tikhonov & Samarskii
  % Input Argument:
  %    - tmesh: Grid of time values, or a vector [t_0, t_final]
  %    - k1, c1, latentHeat: Scalars representing constants in the example. Default: [1, -1, 1]
  %    - tShift: Time shift applied to model problem to avoid
  %      singular data.
  %
  % Output Arguments:
  %    - u_true_T: function u(x, t_final)
  %    - u_true_0: function u(x, t_0)
  %    - u_true_S: function u(s(t), t)
  %    - aux_true_0: function g(t) = a(t) u_x(0, t)
  %    - s_star: s(t_final)
  %    - s_true: function x=s(t)
  %    - u_true: function u=u(x,t)
  %    - a_true: function a=a(t)

  if ~exist('k1', 'var')
    k1 = 1;
  end
  if ~exist('c1', 'var')
    c1 = -1;
  end
  if ~exist('latentHeat', 'var')
    latentHeat = 1;
  end
  if ~exist('tShift', 'var')
    tShift = 1e-1;
  end

  persistent boundary_constant
  persistent B1

  t_initial = tmesh(1);

  t_final = tmesh(end);

  % Analytic diffusion coefficient a(t)
  a_true = @(t) k1*ones(size(t));

  % Equation defining boundary_constant
  if isempty(boundary_constant)
    boundary_constant_f = @(z) (z*latentHeat * sqrt(pi))/2 + (sqrt(k1)*c1*exp(-(z.^2) / (4*k1))) / erf(z / (2 * sqrt(k1)));
    % Assume that the root is in the interval (1e-3, 5)
    boundary_constant = fzero(boundary_constant_f, [1e-3, 5]);
  end

  % Once boundary_constant is known, B1 is fixed.
  if isempty(B1)
    B1 = -c1/erf(boundary_constant/(2*sqrt(k1)));
  end

  % true s(t)
  s_true = @(t) boundary_constant * sqrt(t+tShift);

  % s(T) - s(T) at final time
  s_star = s_true(t_final);

  % Analytic function u=u(x,t)
  u_true = @(x,t) c1 + B1*erf(x./(2*sqrt(k1*(t + tShift))));

  % Analytic function u_x(x,t)
  ux_true = @(x,t) (B1/sqrt(pi*k1*(t + tShift))).*exp(-(x.^2)./(4*k1*(t+tShift)));

  % Analytic function g(t) = a u_x(0,t).
  % This would easily vectorize by writing ux_true(0,t) as a special case
  aux_true_0_ptwise = @(t) a_true(t) .* ux_true(0, t); % == B1*sqrt(k1/(pi*(t + tShift)))
  aux_true_0 = @(t) arrayfun(aux_true_0_ptwise, t);

  % \mu(t) = u(s(t),t) Phase Transition temperature
  u_true_S = @(t) zeros(size(t)); % == u_true(s_true(t), t);

  % w(x) = u(x, T) measurement at final moment
  u_true_T_ptwise = @(x) u_true(x, t_final);
  u_true_T = @(x) arrayfun(u_true_T_ptwise, x);

  % phi(x) = u(x, 0) measurement at initial moment
  u_true_0_ptwise = @(x) u_true(x, t_initial);
  u_true_0 = @(x) arrayfun(u_true_0_ptwise, x);
end
