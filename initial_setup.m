function [max_step_size, u_true_0, mu_meas, w_meas, g, s_star, s_ini, a_ini] = initial_setup(tmesh, xmesh, use_synthetic_data, initial_data_parameter)
% initial_setup: Create parameters and measurements for model problem
% Set use_synthetic_data=true to generate synthetic data for measurements
% Set initial_data_parameter == 0 to have s_ini == s_true
%                            == 1 to have s_ini == s_linear

% The true solution values are used only for the measurements s_star
[u_true_T, u_true_0, u_true_S, g, s_star, s_true, ~, a_true] = true_solution(tmesh);

% Extract initial time value as t_initial
t_initial = tmesh(1);

s0 = s_true(t_initial);

% Extract final time value as t_final
t_final = tmesh(end);

% maximum step size in anti-gradient direction
max_step_size = 1;


% Initial Guess
 %s_linear = @(t) ((s_star - s0)/(t_final - t_initial)) * (t - t_initial) + s0;
 %s_ini = @(t) initial_data_parameter * s_linear(t) + (1 - initial_data_parameter) * s_true(t);
 %s_ini=@(t)((s_star - s0)/(t_final - t_initial)) * (t - t_initial) + s0 + 1/3;
 s_ini=@(t) s_true(t);

% a_ini = @(t) a_true(t);
a_ini = @(t) a_true(t)-0.8;

% Use synthetic data to set measurements
if use_synthetic_data
    svals_ = s_true(tmesh);
    avals_ = a_true(tmesh);
    
    [~, ~ , u_S_, w_meas, ~] = Forward(xmesh, tmesh, svals_, avals_, g, u_true_0);
    
    % Interpolate u_S_ values to form u_true_S
    mu_meas = @(t) interp1(tmesh, u_S_, t);
    disp(['S trace declination: ' num2str(norm(mu_meas(tmesh) - u_true_S(tmesh)))]);
    
    disp(['T trace declination: ' num2str(norm(w_meas(xmesh) - u_true_T(xmesh * svals_(end))))]);
else
    % Use true solution to set measurements
    mu_meas = u_true_S;
    w_meas = u_true_T;
end

end
