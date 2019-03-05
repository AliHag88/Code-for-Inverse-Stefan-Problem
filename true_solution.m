function [u_true_T, u_true_0, u_true_S, g, s_star, s_true, u_true, a_true] = true_solution(tmesh)
  % true_solution: Return parameter functions for "true" solution
  % Input Argument:
  %    - tmesh: Grid of time values, or a vector [t_0, t_final]
  %
  % Output Arguments:
  %    - u_true_T: function u(x, t_final)
  %    - u_true_0: function u(x, t_0)
  %    - u_true_S: function u(s(t), t)
  %    - g: function g(t) = a(t) u_x(0, t)
  %    - s_star: s(t_final)
  %    - s_true: function x=s(t)
  %    - u_true: function u=u(x,t)
  %    - a_true: function a=a(t)

t_initial = tmesh(1);

t_final = tmesh(end);

% Analytic diffusion coefficient a(t)
a_true = @(t) ones(size(t));

% true s(t)
s_true = @(t) (t.^2).*sqrt(t+0.2)+1;

% s(T) - s(T) at final time
s_star = s_true(t_final);

% Analytic function u=u(x,t)
u_true = @(x,t) t .* x .* (1-x);

% Analytic function u_x(x,t)
ux_true = @(x,t) t .* (1-2x);

% Analytic function g(t) = a u_x
ux_true_0 = @(t) a_true(t) .* ux_true(0, t);

% \mu(t) = u(s(t),t) Phase Transition temperature
u_true_S = @(t) u_true(s_true(t), t);

% w(x) = u(x, T) measurement at final moment
u_true_T = @(x) u_true(x, t_final);

% phi(x) = u(x, 0) measurement at initial moment
u_true_0 = @(x) u_true(x, t_initial);
end
