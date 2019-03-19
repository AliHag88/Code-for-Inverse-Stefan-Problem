function [step_size, u_true_0, mu_meas, w_meas, g, s_star, s_ini, a_ini] = initial_setup(tmesh, xmesh)
% The true solution values are used only for the measurements s_star
[u_true_T, u_true_0, u_true_S, g, s_star, s_true, ~, a_true] = true_solution(tmesh);

% Extract initial time value as t_initial
t_initial = tmesh(1);

s0 = s_true(t_initial);

% Extract final time value as t_final
t_final = tmesh(end);

% step_size (size of step in anti-gradient direction)
step_size = 1/300;

% Use true solution to set measurements
%mu_meas = u_true_S;
%w_meas = u_true_T;

% Initial Guess
%s_ini = @(t) (s_star - s0)/(t_final - t_initial) * (t - t_initial) + s0;
s_ini = @(t) s_true(t);
a_ini = a_true;


% Synthetic measurements
svals_ = s_true(tmesh);
avals_ = a_true(tmesh);

[~, ~ , u_S_, u_T_, ~] = Forward(xmesh, tmesh, svals_, avals_, g, u_true_0);

% Interpolate a values to form u_true_S
mu_meas = @(t) interp1(tmesh, u_S_, t);
%mu_meas(tmesh) - u_true_S(tmesh)

% Interpolate w values to form u_true_T
w_meas = @(x) interp1(xmesh, u_T_, x);
%w_meas(xmesh) - u_true_T(xmesh * svals_(end))
end
