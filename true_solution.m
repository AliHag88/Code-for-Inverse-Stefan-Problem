function [u_true_T, u_true_0, u_true_S, s_star, s_true, u_true] = true_solution(tmesh)
  % true_solution: Return parameter functions for "true" solution
  % Input Argument:
  %    - tmesh: Grid of time values, or a vector [t_0, t_final]
  %
  % Output Arguments:
  %    - u_true_T: function u(x, t_final)
  %    - u_true_0: function u(x, t_0)
  %    - u_true_S: function u(s(t), t)
  %    - s_star: s(t_final)
  %    - s_true: function x=s(t)
  %    - u_true: function u=u(x,t)

t_initial = tmesh(1);

t_final = tmesh(end);

% true s(t)
s_true=@(t) (t.^2).*sqrt(t+0.2)+1;

% s(T) - s(T) at final time
s_star = s_true(t_final);

% Analytic function u=u(x,t)
u_true = @(x,t) t .* x .* (1-x);

% \mu(t) = u(s(t),t) Phase Transition temperature
u_true_S = @(t) u_true(s_true(t), t);

% w(x) = u(x, T) measurement at final moment
u_true_T = @(x) u_true(x, t_final);

% phi(x) = u(x, 0) measurement at initial moment
u_true_0 = @(x) u_true(x, t_initial);
end
