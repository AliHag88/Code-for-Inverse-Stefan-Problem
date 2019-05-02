function [max_step_size, u_true_0, mu_meas, w_meas, g, s_star, s_ini, a_ini] = initial_setup(tmesh, xmesh, use_synthetic_data, initial_data_parameter_s, initial_data_parameter_a)
  % initial_setup creates parameters and measurements for model problem
  % Set use_synthetic_data=true to generate synthetic data for measurements
  % Set initial_data_parameter_s == 0 to have s_ini == s_true
  %                              == 1 to have s_ini == s_linear
  % Set initial_data_parameter_a == 0 to have a_ini == a_true.
  % Choosing a nonzero value for this parameter will give a perturbation of the
  % chosen magnitude to a_true.

  % The true solution values are used only for the measurements s_star
  [u_true_T, u_true_0, u_true_S, g, s_star, s_true, ~, a_true] = true_solution(tmesh);

  % Extract initial time value as t_initial
  t_initial = tmesh(1);

  % Calculate (fixed) initial boundary position
  s0 = s_true(t_initial);

  % Extract final time value as t_final
  t_final = tmesh(end);

  % maximum step size in anti-gradient direction
  max_step_size = 1;

  % Initial guess for s(t)
  s_linear = @(t) ((s_star - s0)/(t_final - t_initial)) * (t - t_initial) + s0;
  s_ini = @(t) initial_data_parameter_s * s_linear(t) + (1 - initial_data_parameter_s) * s_true(t);

  % Initial guess for a(t)
  a_ini = @(t) a_true(t) + (4/(t_final - t_initial)) * initial_data_parameter_a * ;

  % Use synthetic data to set measurements
  if use_synthetic_data
    svals_ = s_true(tmesh);
    avals_ = a_true(tmesh);

    [~, ~ , u_S_, w_meas, ~] = Forward(xmesh, tmesh, svals_, avals_, g, u_true_0);

    % Interpolate u_S_ values to form u_true_S
    mu_meas = @(t) interp1(tmesh, u_S_, t);
    fprintf('S trace declination: %2.5f.\n', norm(mu_meas(tmesh) - u_true_S(tmesh)));
    fprintf('T trace declination: %2.5f.\n', norm(w_meas(xmesh) - u_true_T(xmesh * svals_(end))));
  else % Not using synthetic data
    % Use true solution to set measurements
    mu_meas = u_true_S;
    w_meas = u_true_T;
  end % if use_synthetic_data

end % Function
